"""
Common solar physics coordinate systems.

This submodule implements various solar physics coordinate frames for use with
the `astropy.coordinates` module.
"""
import os
import re
import traceback

import numpy as np

import astropy.units as u
from astropy.constants import R_earth
from astropy.coordinates import Attribute, ConvertError, Latitude, Longitude, QuantityAttribute
from astropy.coordinates.baseframe import BaseCoordinateFrame, RepresentationMapping
from astropy.coordinates.representation import (
    CartesianDifferential,
    CartesianRepresentation,
    CylindricalRepresentation,
    SphericalDifferential,
    SphericalRepresentation,
    UnitSphericalRepresentation,
)
from astropy.time import Time
from astropy.utils.data import download_file

from sunpy import log
from sunpy.sun.constants import radius as _RSUN
from sunpy.time.time import _variables_for_parse_time_docstring
from sunpy.util.decorators import add_common_docstring, deprecated, sunpycontextmanager
from sunpy.util.exceptions import warn_user
from .frameattributes import ObserverCoordinateAttribute, TimeFrameAttributeSunPy

_J2000 = Time('J2000.0', scale='tt')

__all__ = ['SunPyBaseCoordinateFrame', 'BaseHeliographic', 'BaseMagnetic',
           'HeliographicStonyhurst', 'HeliographicCarrington',
           'Heliocentric', 'Helioprojective', 'HelioprojectiveRadial',
           'HeliocentricEarthEcliptic', 'GeocentricSolarEcliptic',
           'HeliocentricInertial', 'GeocentricEarthEquatorial',
           'Geomagnetic', 'SolarMagnetic', 'GeocentricSolarMagnetospheric']


def _frame_parameters():
    """
    Returns formatting dictionary to use with add_common_docstring to populate frame docstrings
    """
    ret = {}

    ret['data'] = ("data : `~astropy.coordinates.BaseRepresentation` or ``None``\n"
                   "    A representation object or ``None`` to have no data\n"
                   "    (or use the coordinate component arguments, see below).")
    ret['common'] = (f"obstime : {_variables_for_parse_time_docstring()['parse_time_types']}\n"
                     "    The time of the observation. This is used to determine the\n"
                     "    position of solar-system bodies (e.g., the Sun and the Earth) as\n"
                     "    needed to define the origin and orientation of the frame.\n"
                     "representation_type : `~astropy.coordinates.BaseRepresentation`, str, optional\n"
                     "    A representation class or string name of a representation class.\n"
                     "    This may change the valid coordinate component arguments from the\n"
                     "    defaults (see above). For example, passing\n"
                     "    ``representation_type='cartesian'`` will make the frame expect\n"
                     "    Cartesian coordinate component arguments (typically, ``x``, ``y``,\n"
                     "    and ``z``).\n"
                     "copy : bool, optional\n"
                     "    If `True` (default), make copies of the input coordinate arrays.")
    ret['lonlat'] = ("lon : `~astropy.coordinates.Angle` or `~astropy.units.Quantity`, optional\n"
                     "    The longitude coordinate for this object (``lat`` must also be\n"
                     "    given and ``data`` must be ``None``).\n"
                     "    Not needed if ``data`` is given.\n"
                     "lat : `~astropy.coordinates.Angle` or `~astropy.units.Quantity`, optional\n"
                     "    The latitude coordinate for this object (``lon`` must also be\n"
                     "    given and ``data`` must be ``None``).\n"
                     "    Not needed if ``data`` is given.")
    ret['radius'] = ("radius : `~astropy.units.Quantity`, optional\n"
                     "    The radial distance coordinate from Sun center for this object.\n"
                     "    Defaults to the radius of the Sun. Not needed if ``data`` is given.")
    ret['distance_sun'] = ("distance : `~astropy.units.Quantity`, optional\n"
                           "    The distance coordinate from Sun center for this object.\n"
                           "    Not needed if ``data`` is given.")
    ret['distance_earth'] = ("distance : `~astropy.units.Quantity`, optional\n"
                             "    The distance coordinate from Earth center for this object.\n"
                             "    Not needed if ``data`` is given.")
    ret['xyz'] = ("x : `~astropy.units.Quantity`, optional\n"
                  "    X-axis coordinate for this object. Not needed if ``data`` is given.\n"
                  "y : `~astropy.units.Quantity`, optional\n"
                  "    Y-axis coordinate for this object. Not needed if ``data`` is given.\n"
                  "z : `~astropy.units.Quantity`, optional\n"
                  "    Z-axis coordinate for this object. Not needed if ``data`` is given.")
    ret['observer'] = ("observer : `~sunpy.coordinates.frames.HeliographicStonyhurst`, str\n"
                       "    The location of the observer. If a string is provided,\n"
                       "    it must be a solar system body that can be parsed by\n"
                       "    `~sunpy.coordinates.ephemeris.get_body_heliographic_stonyhurst`\n"
                       "    at the time ``obstime``. Defaults to Earth center.")
    ret['rsun'] = ("rsun : `~astropy.units.Quantity`\n"
                   "    The radius of the Sun in length units. Used to convert a 2D\n"
                   "    coordinate (i.e., no ``radius`` component) to a 3D coordinate by\n"
                   "    assuming that the coordinate is on the surface of the Sun. Defaults\n"
                   "    to the photospheric radius as defined in `sunpy.sun.constants`.")
    ret['equinox'] = (f"equinox : {_variables_for_parse_time_docstring()['parse_time_types']}\n"
                      "    The date for the mean vernal equinox.\n"
                      "    Defaults to the J2000.0 equinox.")
    ret['magnetic_model'] = ("magnetic_model : `str`\n"
                             "    The IGRF model to use for determining the orientation of\n"
                             "    Earth's magnetic dipole pole. The supported options are\n"
                             "    ``'igrf13'`` (default), ``'igrf12'``, ``'igrf11'``, and\n"
                             "    ``'igrf10'``.")
    ret['igrf_reference'] = ("* `International Geomagnetic Reference Field (IGRF) "
                             "<https://www.ngdc.noaa.gov/IAGA/vmod/igrf.html>`__")

    return ret


class SunPyBaseCoordinateFrame(BaseCoordinateFrame):
    """
    Base class for sunpy coordinate frames.

    This class is not intended to be used directly and has no transformations defined.

    * Defines the frame attribute ``obstime`` for observation time.
    * Defines a default wrap angle of 180 degrees for longitude in spherical coordinates,
      which can be overridden via the class variable ``_wrap_angle``.
    * Inject a nice way of representing the object which the coordinate represents.
    """
    obstime = TimeFrameAttributeSunPy()

    default_representation = SphericalRepresentation
    default_differential = SphericalDifferential

    frame_specific_representation_info = {
        SphericalDifferential: [RepresentationMapping('d_lon', 'd_lon', u.arcsec/u.s),
                                RepresentationMapping('d_lat', 'd_lat', u.arcsec/u.s),
                                RepresentationMapping('d_distance', 'd_distance', u.km/u.s)],
    }

    _wrap_angle = 180*u.deg  # for longitude in spherical coordinates

    def __init__(self, *args, **kwargs):
        self.object_name = None

        # If wrap_longitude=False is passed in, do not impose a specific wrap angle for the frame
        if not kwargs.pop('wrap_longitude', True):
            self._wrap_angle = None

        # If obstime is not provided but observer has an obstime, use that as the obstime
        if 'obstime' not in kwargs and 'observer' in kwargs and getattr(kwargs['observer'], 'obstime', None) is not None:
            kwargs['obstime'] = kwargs['observer'].obstime

        super().__init__(*args, **kwargs)

        # If obstime is specified, treat the default observer (None) as explicitly set
        if self.obstime is not None and self.is_frame_attr_default('observer'):
            self._attr_names_with_defaults.remove('observer')

        return

    def represent_as(self, base, s='base', in_frame_units=False):
        data = super().represent_as(base, s, in_frame_units=in_frame_units)

        # If a frame wrap angle is set, use that wrap angle for any spherical representations.
        if self._wrap_angle is not None and \
           isinstance(data, UnitSphericalRepresentation | SphericalRepresentation):
            data.lon.wrap_angle = self._wrap_angle
        return data

    def __str__(self):
        # We override this here so that when you print a SkyCoord it shows the
        # observer as the string and not the whole massive coordinate.
        if getattr(self, "object_name", None):
            return f"<{self.__class__.__name__} Coordinate for '{self.object_name}'>"
        else:
            return super().__str__()

    @property
    def _is_2d(self):
        return (self._data is not None and self._data.norm().unit is u.one
                and u.allclose(self._data.norm(), 1*u.one))


class BaseHeliographic(SunPyBaseCoordinateFrame):
    """
    Base class for HeliographicCarrington (HGC) and HeliographicStonyhurst (HGS) frames.

    This class is not intended to be used directly and has no transformations defined.
    """
    frame_specific_representation_info = {
        SphericalRepresentation: [RepresentationMapping('lon', 'lon', u.deg),
                                  RepresentationMapping('lat', 'lat', u.deg),
                                  RepresentationMapping('distance', 'radius', None)],
        SphericalDifferential: [RepresentationMapping('d_lon', 'd_lon', u.arcsec/u.s),
                                RepresentationMapping('d_lat', 'd_lat', u.arcsec/u.s),
                                RepresentationMapping('d_distance', 'd_radius', u.km/u.s)],
    }

    rsun = QuantityAttribute(default=_RSUN, unit=u.km)

    def make_3d(self):
        """
        Returns a fully 3D coordinate based on this coordinate.

        If this coordinate is only 2D (i.e., no ``radius`` component) or is a
        unit vector (i.e., the norm of the coordinate is unity), a new
        coordinate is created that corresponds to the surface of the Sun.
        That is, the 3D coordinate will retain the ``lon`` and ``lat``, and
        ``radius`` will be set to the frame's ``rsun`` frame attribute.

        If this coordinate is already fully 3D, it is directly returned, even
        if it does not lie on the surface of the Sun.

        Returns
        -------
        frame : `~sunpy.coordinates.frames.BaseHeliographic`
            The fully 3D coordinate
        """
        if self._is_2d:
            return self.realize_frame(self._data * self.rsun)

        # The coordinate is already 3D
        return self


@add_common_docstring(**_frame_parameters())
class HeliographicStonyhurst(BaseHeliographic):
    """
    A coordinate or frame in the Stonyhurst Heliographic (HGS) system.

    - The origin is the center of the Sun.
    - The Z-axis (+90 degrees latitude) is aligned with the Sun's north pole.
    - The X-axis (0 degrees longitude and 0 degrees latitude) is aligned with the projection of
      the Sun-Earth line onto the Sun's equatorial plane.

    This system is also know as the Heliocentric Earth Equatorial (HEEQ) system when
    represented using Cartesian components.

    A new instance can be created using the following signatures
    (note that if supplied, ``obstime`` and ``representation_type`` must be
    keyword arguments)::

        HeliographicStonyhurst(lon, lat, obstime=obstime)
        HeliographicStonyhurst(lon, lat, radius, obstime=obstime)
        HeliographicStonyhurst(x, y, z, representation_type='cartesian', obstime=obstime)

    Parameters
    ----------
    {data}
    {lonlat}
    {radius}
    {rsun}
    {common}

    Examples
    --------
    >>> from astropy.coordinates import SkyCoord
    >>> import sunpy.coordinates
    >>> import astropy.units as u
    >>> sc = SkyCoord(1*u.deg, 1*u.deg, 2*u.km,
    ...               frame="heliographic_stonyhurst",
    ...               obstime="2010/01/01T00:00:45")
    >>> sc
    <SkyCoord (HeliographicStonyhurst: obstime=2010-01-01T00:00:45.000, rsun=695700.0 km): (lon, lat, radius) in (deg, deg, km)
        (1., 1., 2.)>
    >>> sc.frame
    <HeliographicStonyhurst Coordinate (obstime=2010-01-01T00:00:45.000, rsun=695700.0 km): (lon, lat, radius) in (deg, deg, km)
        (1., 1., 2.)>
    >>> sc = SkyCoord(HeliographicStonyhurst(-10*u.deg, 2*u.deg))
    >>> sc
    <SkyCoord (HeliographicStonyhurst: obstime=None, rsun=695700.0 km): (lon, lat) in deg
        (-10., 2.)>
    >>> sc = SkyCoord(CartesianRepresentation(0*u.km, 45*u.km, 2*u.km),
    ...               obstime="2011/01/05T00:00:50",
    ...               frame="heliographic_stonyhurst")
    >>> sc
    <SkyCoord (HeliographicStonyhurst: obstime=2011-01-05T00:00:50.000, rsun=695700.0 km): (lon, lat, radius) in (deg, deg, km)
    (90., 2.54480438, 45.04442252)>
    """
    name = "heliographic_stonyhurst"

    def _apply_diffrot(self, duration, rotation_model):
        oldrepr = self.spherical

        from sunpy.sun.models import differential_rotation
        log.debug(f"Applying {duration} of solar rotation")
        newlon = oldrepr.lon + differential_rotation(duration,
                                                     oldrepr.lat,
                                                     model=rotation_model,
                                                     frame_time='sidereal')
        newrepr = SphericalRepresentation(newlon, oldrepr.lat, oldrepr.distance)

        return self.realize_frame(newrepr)


@add_common_docstring(**_frame_parameters())
class HeliographicCarrington(BaseHeliographic):
    """
    A coordinate or frame in the Carrington Heliographic (HGC) system.

    - The origin is the center of the Sun.
    - The Z-axis (+90 degrees latitude) is aligned with the Sun's north pole.
    - The X-axis and Y-axis rotate with a period of 25.38 days.

    This system differs from Stonyhurst Heliographic (HGS) in its definition of longitude. This
    longitude is an "apparent" longitude because it takes into account the time it takes for light
    to travel from the Sun's surface to the observer (see :ref:`sunpy-topic-guide-coordinates-carrington`).
    Thus, the observer needs to be specified to be able to transform to any other coordinate frame.

    A new instance can be created using the following signatures
    (note that if supplied, ``obstime`` and ``observer`` must be a keyword argument)::

        HeliographicCarrington(lon, lat, obstime=obstime, observer=observer)
        HeliographicCarrington(lon, lat, radius, obstime=obstime, observer=observer)

    If you want to define the location in HGC such that the observer for the coordinate frame is
    the same as that location (e.g., the location of an observatory in its corresponding HGC
    frame), use ``observer='self'``::

        HeliographicCarrington(lon, lat, radius, obstime=obstime, observer='self')

    Parameters
    ----------
    {data}
    {lonlat}
    {radius}
    {observer}
    {rsun}
    {common}

    Examples
    --------
    >>> from astropy.coordinates import SkyCoord
    >>> import sunpy.coordinates
    >>> import astropy.units as u
    >>> sc = SkyCoord(1*u.deg, 2*u.deg, 3*u.km,
    ...               frame="heliographic_carrington",
    ...               observer="earth",
    ...               obstime="2010/01/01T00:00:30")
    >>> sc
    <SkyCoord (HeliographicCarrington: obstime=2010-01-01T00:00:30.000, rsun=695700.0 km, observer=<HeliographicStonyhurst Coordinate for 'earth'>): (lon, lat, radius) in (deg, deg, km)
        (1., 2., 3.)>

    >>> sc = SkyCoord([1,2,3]*u.deg, [4,5,6]*u.deg, [5,6,7]*u.km,
    ...               obstime="2010/01/01T00:00:45",
    ...               observer="self",
    ...               frame="heliographic_carrington")
    >>> sc
    <SkyCoord (HeliographicCarrington: obstime=2010-01-01T00:00:45.000, rsun=695700.0 km, observer=self): (lon, lat, radius) in (deg, deg, km)
        [(1., 4., 5.), (2., 5., 6.), (3., 6., 7.)]>

    >>> sc = SkyCoord(CartesianRepresentation(0*u.km, 45*u.km, 2*u.km),
    ...               obstime="2011/01/05T00:00:50",
    ...               frame="heliographic_carrington")
    >>> sc
    <SkyCoord (HeliographicCarrington: obstime=2011-01-05T00:00:50.000, rsun=695700.0 km, observer=None): (lon, lat, radius) in (deg, deg, km)
        (90., 2.54480438, 45.04442252)>
    """
    name = "heliographic_carrington"
    _wrap_angle = 360*u.deg

    observer = ObserverCoordinateAttribute(HeliographicStonyhurst)

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        if not isinstance(self.observer, BaseCoordinateFrame) and self.observer == 'self' and self._is_2d:
            raise ValueError("Full 3D coordinate (including radius) must be specified "
                             "when observer='self'.")


@add_common_docstring(**_frame_parameters())
class Heliocentric(SunPyBaseCoordinateFrame):
    """
    A coordinate or frame in the Heliocentric system, which is observer-based.

    - The origin is the center of the Sun.
    - The Z-axis is aligned with the Sun-observer line.
    - The Y-axis is aligned with the component of the vector to the Sun's north pole that is
      perpendicular to the Z-axis.

    This frame defaults to a Cartesian component representation, which is known as Heliocentric
    Cartesian (HCC). This frame can also be represented using cylindrical components, where
    where ``rho`` is the impact parameter and ``psi`` is the position angle.
    ``psi`` is measured relative to the west limb, rather than solar north, so is shifted
    by 90 degrees compared to the convention of the Heliocentric Radial (HCR) system.

    A new instance can be created using the following signatures
    (note that if supplied, ``obstime``, ``observer``, and ``representation_type`` must be
    keyword arguments)::

        Heliocentric(x, y, z, obstime=obstime, observer=observer)
        Heliocentric(rho, psi, z, representation_type='cylindrical', obstime=obstime, observer=observer)

    Parameters
    ----------
    {data}
    {xyz}
    {observer}
    {common}

    Examples
    --------

    >>> from astropy.coordinates import SkyCoord, CartesianRepresentation
    >>> import sunpy.coordinates
    >>> import astropy.units as u

    >>> sc = SkyCoord(CartesianRepresentation(10*u.km, 1*u.km, 2*u.km),
    ...               obstime="2011/01/05T00:00:50", observer="earth", frame="heliocentric")
    >>> sc
    <SkyCoord (Heliocentric: obstime=2011-01-05T00:00:50.000, observer=<HeliographicStonyhurst Coordinate for 'earth'>): (x, y, z) in km
        (10., 1., 2.)>

    >>> sc = SkyCoord([1,2]*u.km, [3,4]*u.m, [5,6]*u.cm,
    ...               obstime="2011/01/01T00:00:54", observer="earth", frame="heliocentric")
    >>> sc
    <SkyCoord (Heliocentric: obstime=2011-01-01T00:00:54.000, observer=<HeliographicStonyhurst Coordinate for 'earth'>): (x, y, z) in (km, m, cm)
        [(1., 3., 5.), (2., 4., 6.)]>

    >>> sc = SkyCoord(CylindricalRepresentation(10*u.km, 60*u.deg, 10*u.km),
    ...               obstime="2011/01/05T00:00:50", observer="earth", frame="heliocentric")
    >>> sc
    <SkyCoord (Heliocentric: obstime=2011-01-05T00:00:50.000, observer=<HeliographicStonyhurst Coordinate for 'earth'>): (x, y, z) in km
        (5., 8.66025404, 10.)>
    """
    default_representation = CartesianRepresentation
    default_differential = CartesianDifferential

    frame_specific_representation_info = {
        CylindricalRepresentation: [RepresentationMapping('phi', 'psi', u.deg)]
    }

    observer = ObserverCoordinateAttribute(HeliographicStonyhurst)

    def represent_as(self, base, s='base', in_frame_units=False):
        data = super().represent_as(base, s, in_frame_units=in_frame_units)

        # For cylindrical representations, wrap the `psi` component (natively `phi`) at 360 deg
        if isinstance(data, CylindricalRepresentation):
            data.phi.wrap_at(360*u.deg, inplace=True)
        return data


@add_common_docstring(**_frame_parameters())
class Helioprojective(SunPyBaseCoordinateFrame):
    """
    A coordinate or frame in the Helioprojective Cartesian (HPC) system.

    This is an observer-based spherical coordinate system, with:

    - The origin is the observer location.
    - The line connecting the poles of the frame is parallel to the component of the Sun's
      rotation axis that is perpendicular to the observer-Sun line.
    - ``Tx`` (short for "theta_x", or :math:`\\theta_x`) is the longitude, the angle relative to
      the plane containing the observer-Sun line and the poles of the frame, with positive values
      in the direction of the Sun's west limb.
    - ``Ty`` (short for "theta_y", or :math:`\\theta_y`) is the latitude, the angle relative to the
      plane that is perpendicular to the poles of the frame, with positive values in the direction
      of the Sun's north pole.
    - ``distance`` is the observer-object distance.

    .. note::
        It can be confusing that the name of this coordinate system uses the adjective "Cartesian"
        despite being a spherical coordinate system. The reason for this name is because close to
        the center of the Sun where the small-angle approximation is appropriate, ``Tx`` and ``Ty``
        can be thought of as two components of a Cartesian-like coordinate system. The corresponding
        third Cartesian-like component, sometimes called zeta (:math:`\\zeta`), is defined to be the
        Sun-observer distance (``observer.radius``) minus the observer-object distance (``distance``).

    This system is frequently used in a projective form without ``distance`` specified. When
    transforming such a 2D coordinate to another frame, the object location usually needs to be
    fully 3D, which is achieved by calling :meth:`.make_3d` to generate the ``distance`` component.
    The default assumption is that the object lies on the surface of the Sun if the 2D coordinate is
    on the solar disk, but is otherwise undefined if the 2D coordinate is beyond the solar limb.
    This assumption can be modified using a screen (e.g., :func:`~sunpy.coordinates.SphericalScreen`).

    A new instance can be created using the following signatures
    (note that if supplied, ``obstime`` and ``observer`` must be keyword arguments)::

        Helioprojective(Tx, Ty, obstime=obstime, observer=observer)
        Helioprojective(Tx, Ty, distance, obstime=obstime, observer=observer)

    Parameters
    ----------
    {data}
    Tx : `~astropy.coordinates.Angle` or `~astropy.units.Quantity`
        The theta_x (:math:`\\theta_x`) component for this object. Not needed if ``data`` is given.
    Ty : `~astropy.coordinates.Angle` or `~astropy.units.Quantity`
        The theta_y (:math:`\\theta_y`) component for this object. Not needed if ``data`` is given.
    distance : `~astropy.units.Quantity`
        The distance component from the observer for this object. Not needed if ``data`` is given.
    {observer}
    {rsun}
    {common}

    See Also
    --------
    HelioprojectiveRadial
    ~sunpy.coordinates.PlanarScreen, ~sunpy.coordinates.SphericalScreen

    Examples
    --------
    >>> from astropy.coordinates import SkyCoord
    >>> import sunpy.coordinates
    >>> import astropy.units as u
    >>> sc = SkyCoord(0*u.deg, 0*u.deg, 5*u.km,
    ...               obstime="2010/01/01T00:00:00", observer="earth", frame="helioprojective")
    >>> sc
    <SkyCoord (Helioprojective: obstime=2010-01-01T00:00:00.000, rsun=695700.0 km, observer=<HeliographicStonyhurst Coordinate for 'earth'>): (Tx, Ty, distance) in (arcsec, arcsec, km)
        (0., 0., 5.)>
    >>> sc = SkyCoord(0*u.deg, 0*u.deg,
    ...               obstime="2010/01/01T00:00:00", observer="earth", frame="helioprojective")
    >>> sc
    <SkyCoord (Helioprojective: obstime=2010-01-01T00:00:00.000, rsun=695700.0 km, observer=<HeliographicStonyhurst Coordinate for 'earth'>): (Tx, Ty) in arcsec
        (0., 0.)>
    >>> sc = SkyCoord(CartesianRepresentation(1*u.AU, 1e5*u.km, -2e5*u.km),
    ...               obstime="2011/01/05T00:00:50", observer="earth", frame="helioprojective")
    >>> sc
    <SkyCoord (Helioprojective: obstime=2011-01-05T00:00:50.000, rsun=695700.0 km, observer=<HeliographicStonyhurst Coordinate for 'earth'>): (Tx, Ty, distance) in (arcsec, arcsec, AU)
        (137.87948623, -275.75878762, 1.00000112)>
    """
    frame_specific_representation_info = {
        SphericalRepresentation: [RepresentationMapping('lon', 'Tx', u.arcsec),
                                  RepresentationMapping('lat', 'Ty', u.arcsec),
                                  RepresentationMapping('distance', 'distance', None)],
        SphericalDifferential: [RepresentationMapping('d_lon', 'd_Tx', u.arcsec/u.s),
                                RepresentationMapping('d_lat', 'd_Ty', u.arcsec/u.s),
                                RepresentationMapping('d_distance', 'd_distance', u.km/u.s)],
        UnitSphericalRepresentation: [RepresentationMapping('lon', 'Tx', u.arcsec),
                                      RepresentationMapping('lat', 'Ty', u.arcsec)],
    }

    rsun = QuantityAttribute(default=_RSUN, unit=u.km)
    observer = ObserverCoordinateAttribute(HeliographicStonyhurst)

    @property
    def angular_radius(self):
        """
        Angular radius of the Sun as seen by the observer.

        The ``rsun`` frame attribute is the radius of the Sun in length units.
        The tangent vector from the observer to the edge of the Sun forms a
        right-angle triangle with the radius of the Sun as the far side and the
        Sun-observer distance as the hypotenuse. Thus, the sine of the angular
        radius of the Sun is ratio of these two distances.
        """
        from sunpy.coordinates.sun import _angular_radius  # avoiding a circular import

        if not isinstance(self.observer, HeliographicStonyhurst):
            if self.observer is None:
                raise ValueError("The observer must be defined, not `None`.")
            raise ValueError("The observer must be fully defined by specifying `obstime`.")
        return _angular_radius(self.rsun, self.observer.radius)

    def make_3d(self):
        """
        This method calculates the third coordinate of the Helioprojective
        frame. It assumes that the coordinate point is on the surface of the Sun.

        If a point in the frame is off limb then NaN will be returned.

        Returns
        -------
        new_frame : `~sunpy.coordinates.frames.Helioprojective`
            A new frame instance with all the attributes of the original but
            now with a third coordinate.
        """
        # Skip if we already are 3D
        if not self._is_2d:
            return self

        if not isinstance(self.observer, BaseCoordinateFrame):
            raise ConvertError("Cannot calculate distance to the Sun "
                               f"for observer '{self.observer}' "
                               "without `obstime` being specified.")

        rep = self.represent_as(UnitSphericalRepresentation)
        lat, lon = rep.lat, rep.lon

        # Check for the use of floats with lower precision than the native Python float
        if not set([lon.dtype.type, lat.dtype.type]).issubset([float, np.float64, np.longdouble]):
            warn_user("The Helioprojective component values appear to be lower "
                      "precision than the native Python float: "
                      f"Tx is {lon.dtype.name}, and Ty is {lat.dtype.name}. "
                      "To minimize precision loss, you may want to cast the values to "
                      "`float` or `numpy.float64` via the NumPy method `.astype()`.")

        # Calculate the distance to the surface of the Sun using the law of cosines
        cos_alpha = np.cos(lat) * np.cos(lon)
        c = self.observer.radius**2 - self.rsun**2
        b = -2 * self.observer.radius * cos_alpha
        # Ignore sqrt of NaNs
        with np.errstate(invalid='ignore'):
            d = ((-1*b) - np.sqrt(b**2 - 4*c)) / 2  # use the "near" solution

        if self._assumed_screen:
            d_screen = self._assumed_screen.calculate_distance(self)
            d = np.fmin(d, d_screen) if self._assumed_screen.only_off_disk else d_screen

        # This warning can be triggered in specific draw calls when plt.show() is called
        # we can not easily prevent this, so we check the specific function is being called
        # within the stack trace.
        stack_trace = traceback.format_stack()
        matching_string = 'wcsaxes.*(_draw_grid|_update_ticks)'
        bypass = any([re.search(matching_string, string) for string in stack_trace])
        if not bypass and np.all(np.isnan(d)) and np.any(np.isfinite(cos_alpha)):
            warn_user("The conversion of these 2D helioprojective coordinates to 3D is all NaNs "
                      "because off-disk coordinates need an additional assumption to be mapped to "
                      "calculate distance from the observer. Consider using the context manager "
                      "`SphericalScreen()`.")

        return self.realize_frame(SphericalRepresentation(lon=lon,
                                                          lat=lat,
                                                          distance=d))

    @u.quantity_input
    def is_visible(self, *, tolerance: u.m = 1*u.m):
        """
        Returns whether the coordinate is on the visible side of the Sun.

        A coordinate is visible if it can been seen from the observer (the ``observer``
        frame attribute) assuming that the Sun is an opaque sphere with a fixed radius
        (the ``rsun`` frame attribute). The visible side of the Sun is always smaller
        than a full hemisphere.

        Parameters
        ----------
        tolerance : `~astropy.units.Quantity`
            The depth below the surface of the Sun that should be treated as
            transparent.

        Notes
        -----
        If the coordinate is 2D, it is automatically deemed visible. A 2D coordinate
        describes a look direction from the observer, who would simply see whatever is
        in "front", and thus cannot correspond to a point hidden from the observer.

        The ``tolerance`` parameter accommodates situations where the limitations of
        numerical precision would falsely conclude that a coordinate is not visible.
        For example, a coordinate that is expressly created to be on the solar surface
        may be calculated to be slightly below the surface, and hence not visible if
        there is no tolerance. However, a consequence of the ``tolerance`` parameter
        is that a coordinate that is formally on the far side of the Sun but is
        extremely close to the solar limb can be evaluated as visible. With the
        default ``tolerance`` value of 1 meter, a coordinate on the surface of the Sun
        can be up to 11 arcseconds of heliographic longitude past the solar limb and
        still be evaluated as visible.

        Examples
        --------
        >>> import numpy as np
        >>> import astropy.units as u
        >>> from astropy.coordinates import SkyCoord
        >>> from sunpy.coordinates import HeliographicStonyhurst, Helioprojective

        >>> hpc_frame = Helioprojective(observer='earth', obstime='2023-08-03')

        >>> in_front = SkyCoord(0*u.arcsec, 0*u.arcsec, 0.5*u.AU, frame=hpc_frame)
        >>> print(in_front.is_visible())
        True

        >>> behind = SkyCoord(0*u.arcsec, 0*u.arcsec, 1.5*u.AU, frame=hpc_frame)
        >>> print(behind.is_visible())
        False

        >>> hgs_array = SkyCoord(np.arange(-180, 180, 60)*u.deg, [0]*6*u.deg,
        ...                      frame='heliographic_stonyhurst', obstime='2023-08-03')
        >>> print(hgs_array)
        <SkyCoord (HeliographicStonyhurst: obstime=2023-08-03T00:00:00.000, rsun=695700.0 km): (lon, lat) in deg
            [(-180., 0.), (-120., 0.), ( -60., 0.), (   0., 0.), (  60., 0.),
             ( 120., 0.)]>
        >>> print(hgs_array.transform_to(hpc_frame).is_visible())
        [False False  True  True  True False]
        """
        # If the coordinate is 2D, it must be visible
        if self._is_2d:
            return np.ones_like(self.data, dtype=bool)

        # Use a slightly smaller solar radius to accommodate numerical precision
        solar_radius = self.rsun - tolerance

        data = self.cartesian
        data_to_sun = self.observer.radius * CartesianRepresentation(1, 0, 0) - data

        # When representing the helioprojective point as true Cartesian, the X value is the
        # distance from the observer to the point in the sunward direction
        is_behind_observer = data.x < 0
        # When comparing heliocentric angles, we compare the sine values and hence avoid calling arcsin()
        is_beyond_limb = np.sqrt(data.y **2 + data.z **2) / data.norm() > solar_radius / self.observer.radius
        is_above_surface = data_to_sun.norm() >= solar_radius
        is_on_near_side = data.dot(data_to_sun) >= 0

        return is_behind_observer | is_beyond_limb | (is_on_near_side & is_above_surface)

    _assumed_screen = None

    @classmethod
    @sunpycontextmanager
    @deprecated('6.0', alternative='sunpy.coordinates.screens.SphericalScreen')
    def assume_spherical_screen(cls, center, only_off_disk=False, *, radius=None):
        try:
            old_assumed_screen = cls._assumed_screen  # nominally None
            from sunpy.coordinates import SphericalScreen
            sph_screen = SphericalScreen(center, radius=radius, only_off_disk=only_off_disk)
            cls._assumed_screen = sph_screen
            yield
        finally:
            cls._assumed_screen = old_assumed_screen


@add_common_docstring(**_frame_parameters())
class HelioprojectiveRadial(SunPyBaseCoordinateFrame):
    """
    A coordinate or frame in the Helioprojective Radial system.

    This is an observer-based spherical coordinate system, with:

    - ``psi`` is the position angle of the coordinate, measured eastward from solar
      north
    - ``delta`` is the declination angle, which is the impact angle (the angle
      between the observer-Sun line and the observer-coordinate line) minus 90
      degrees
    - ``r`` is the observer-coordinate distance

    .. note::
        The declination angle, rather than the impact angle, is used as a component
        in order to match the FITS WCS definition. The impact angle can be readily
        retrieved using the `theta` property.

    Parameters
    ----------
    {data}
    psi : `~astropy.coordinates.Angle` or `~astropy.units.Quantity`
        The position angle. Not needed if ``data`` is given.
    delta : `~astropy.coordinates.Angle` or `~astropy.units.Quantity`
        The declination angle. Not needed if ``data`` is given.
    r: `~astropy.coordinates.Angle` or `~astropy.units.Quantity`
        The observer-coordinate distance. Not needed if ``data`` is given.
    {rsun}
    {observer}
    {common}

    See Also
    --------
    Helioprojective

    Examples
    --------
    >>> from astropy.coordinates import SkyCoord
    >>> import sunpy.coordinates
    >>> import astropy.units as u

    >>> sc = SkyCoord(0*u.deg, -90*u.deg, 5*u.km,
    ...               obstime="2010/01/01T00:00:00", observer="earth", frame="helioprojectiveradial")
    >>> sc
    <SkyCoord (HelioprojectiveRadial: obstime=2010-01-01T00:00:00.000, rsun=695700.0 km, observer=<HeliographicStonyhurst Coordinate for 'earth'>): (psi, delta, distance) in (deg, deg, km)
        (0., -90., 5.)>
    >>> sc.theta
    <Angle 0. arcsec>

    >>> sc = SkyCoord(30*u.deg, -89.9*u.deg,
    ...               obstime="2010/01/01T00:00:00", observer="earth", frame="helioprojectiveradial")
    >>> sc
    <SkyCoord (HelioprojectiveRadial: obstime=2010-01-01T00:00:00.000, rsun=695700.0 km, observer=<HeliographicStonyhurst Coordinate for 'earth'>): (psi, delta) in deg
        (30., -89.9)>
    >>> sc.theta
    <Angle 360. arcsec>

    >>> sc = SkyCoord(CartesianRepresentation(1e5*u.km, -2e5*u.km, -1*u.AU),
    ...               obstime="2011/01/05T00:00:50", observer="earth", frame="helioprojectiveradial")
    >>> sc
    <SkyCoord (HelioprojectiveRadial: obstime=2011-01-05T00:00:50.000, rsun=695700.0 km, observer=<HeliographicStonyhurst Coordinate for 'earth'>): (psi, delta, distance) in (deg, deg, km)
        (296.56505118, -89.91435897, 1.49598038e+08)>
    >>> sc.theta
    <Angle 308.30772022 arcsec>

    .. minigallery:: sunpy.coordinates.HelioprojectiveRadial
    """
    _wrap_angle = 360*u.deg

    default_representation = SphericalRepresentation

    frame_specific_representation_info = {
        SphericalRepresentation: [RepresentationMapping('lon', 'psi', u.deg),
                                  RepresentationMapping('lat', 'delta', u.deg),
                                  RepresentationMapping('distance', 'distance', None)],
        SphericalDifferential: [RepresentationMapping('d_lon', 'd_psi', u.deg/u.s),
                                RepresentationMapping('d_lat', 'd_delta', u.deg/u.s),
                                RepresentationMapping('d_distance', 'd_distance', u.km/u.s)],
        UnitSphericalRepresentation: [RepresentationMapping('lon', 'psi', u.deg),
                                      RepresentationMapping('lat', 'delta', u.deg)],
    }

    rsun = QuantityAttribute(default=_RSUN, unit=u.km)
    observer = ObserverCoordinateAttribute(HeliographicStonyhurst)

    @property
    def theta(self):
        """
        Returns the impact angle, which is the declination angle plus 90 degrees.
        """
        return (90*u.deg + self.spherical.lat).to(u.arcsec)

    def make_3d(self):
        """
        Returns a 3D version of this coordinate.

        If the coordinate is 2D, the default assumption is that the coordinate is on
        the surface of the Sun, and the distance component is calculated
        accordingly. Under this assumption, if the 2D coordinate is outside the
        disk, the distance component will be NaN.

        The assumption can be changed using one of the screens in
        `sunpy.coordinates.screens`.

        Returns
        -------
        `~sunpy.coordinates.frames.HelioprojectiveRadial`
            The 3D version of this coordinate.
        """
        # Skip if we already are 3D
        if not self._is_2d:
            return self

        # Make 3D by going through HPC, which thus will make use of any screen
        hpc_frame = Helioprojective(obstime=self.obstime, observer=self.observer, rsun=self.rsun)
        return self.transform_to(hpc_frame).make_3d().transform_to(self)


@add_common_docstring(**_frame_parameters())
class HeliocentricEarthEcliptic(SunPyBaseCoordinateFrame):
    """
    A coordinate or frame in the Heliocentric Earth Ecliptic (HEE) system.

    - The origin is the center of the Sun.
    - The X-axis (0 degrees longitude and 0 degrees latitude) is aligned with the Sun-Earth line.
    - The Z-axis (+90 degrees latitude) is aligned with the component perpendicular to the X-axis
      of the mean ecliptic pole at the observation time.

    Parameters
    ----------
    {data}
    {lonlat}
    {distance_sun}
    {common}
    """


@add_common_docstring(**_frame_parameters())
class GeocentricSolarEcliptic(SunPyBaseCoordinateFrame):
    """
    A coordinate or frame in the Geocentric Solar Ecliptic (GSE) system.

    - The origin is the center of the Earth.
    - The X-axis (0 degrees longitude and 0 degrees latitude) is aligned with the Earth-Sun line.
    - The Z-axis (+90 degrees latitude) is aligned with the component perpendicular to the X-axis
      of the mean ecliptic pole at the observation time.

    Parameters
    ----------
    {data}
    {lonlat}
    {distance_earth}
    {common}

    Notes
    -----
    Aberration due to Earth motion is not included.
    """


@add_common_docstring(**_frame_parameters())
class HeliocentricInertial(SunPyBaseCoordinateFrame):
    """
    A coordinate or frame in the Heliocentric Inertial (HCI) system.

    - The origin is the center of the Sun.
    - The Z-axis (+90 degrees latitude) is aligned with the Sun's north pole.
    - The X-axis (0 degrees longitude and 0 degrees latitude) is aligned with the solar ascending
      node on the ecliptic (mean J2000.0).

    Parameters
    ----------
    {data}
    {lonlat}
    {distance_sun}
    {common}

    Notes
    -----
    The solar ascending node on the ecliptic lies on the intersection of the solar equatorial
    plane with the ecliptic plane, not on the intersection of the celestial equatorial plane with
    the ecliptic plane.
    """


@add_common_docstring(**_frame_parameters())
class GeocentricEarthEquatorial(SunPyBaseCoordinateFrame):
    """
    A coordinate or frame in the Geocentric Earth Equatorial (GEI) system.

    - The origin is the center of the Earth.
    - The Z-axis (+90 degrees latitude) is aligned with the Earth's north pole.
    - The X-axis (0 degrees longitude and 0 degrees latitude) is aligned with the mean (not true)
      vernal equinox.

    Parameters
    ----------
    {data}
    {lonlat}
    {distance_earth}
    {equinox}
    {common}

    Notes
    -----
    Aberration due to Earth motion is not included.
    """
    equinox = TimeFrameAttributeSunPy(default=_J2000)


class BaseMagnetic(SunPyBaseCoordinateFrame):
    """
    Base class for frames that rely on the Earth's magnetic model (MAG, SM, and GSM).

    This class is not intended to be used directly and has no transformations defined.
    """
    magnetic_model = Attribute(default='igrf13')

    @property
    def _igrf_file(self):
        if not self.magnetic_model.startswith("igrf"):
            raise ValueError

        # First look if the file is bundled in package
        local_file = os.path.join(os.path.dirname(__file__), "data",
                                  f"{self.magnetic_model}coeffs.txt")
        if os.path.exists(local_file):
            return local_file

        # Otherwise download the file and cache it
        return download_file("https://www.ngdc.noaa.gov/IAGA/vmod/coeffs/"
                             f"{self.magnetic_model}coeffs.txt", cache=True)

    @property
    def _lowest_igrf_coeffs(self):
        with open(self._igrf_file) as f:
            while not (line := f.readline()).startswith('g/h'):
                pass

            years = list(map(float, line.split()[3:-1]))
            g10s = list(map(float, f.readline().split()[3:]))
            g11s = list(map(float, f.readline().split()[3:]))
            h11s = list(map(float, f.readline().split()[3:]))

        decimalyear = self.obstime.utc.decimalyear
        output_shape = decimalyear.shape
        if np.any(decimalyear < 1900.0):
            raise ValueError("At least one of the dates is earlier than the year 1900, which is unsupported.")

        # Use piecewise linear interpolation before the last year
        decimalyear = np.atleast_1d(decimalyear)
        g10 = np.interp(decimalyear, years, g10s[:-1])
        g11 = np.interp(decimalyear, years, g11s[:-1])
        h11 = np.interp(decimalyear, years, h11s[:-1])

        # Use secular variation beyond the last year
        future = decimalyear > years[-1]
        if np.sum(future) > 0:
            g10[future] = g10s[-2] + (decimalyear[future] - years[-1]) * g10s[-1]
            g11[future] = g11s[-2] + (decimalyear[future] - years[-1]) * g11s[-1]
            h11[future] = h11s[-2] + (decimalyear[future] - years[-1]) * h11s[-1]

        g10 = g10.reshape(output_shape)
        g11 = g11.reshape(output_shape)
        h11 = h11.reshape(output_shape)

        return g10, g11, h11

    @add_common_docstring(**_frame_parameters())
    @property
    def dipole_lonlat(self):
        """
        The geographic longitude/latitude of the Earth's magnetic north pole.

        This position is calculated from the first three coefficients of the selected
        IGRF model per Franz & Harper (2002). The small offset between dipole center
        and Earth center is ignored.

        References
        ----------
        {igrf_reference}
        """
        g10, g11, h11 = self._lowest_igrf_coeffs
        # Intentionally use arctan() instead of arctan2() to get angles in specific quadrants
        lon = (np.arctan(h11 / g11) << u.rad).to(u.deg)
        lat = 90*u.deg - np.arctan((g11 * np.cos(lon) + h11 * np.sin(lon)) / g10)
        return Longitude(lon, wrap_angle=180*u.deg), Latitude(lat)

    @add_common_docstring(**_frame_parameters())
    @property
    def dipole_moment(self):
        """
        The Earth's dipole moment.

        The moment is calculated from the first three coefficients of the selected
        IGRF model per Franz & Harper (2002).

        References
        ----------
        {igrf_reference}
        """
        g10, g11, h11 = self._lowest_igrf_coeffs
        moment = np.sqrt(g10**2 + g11**2 + h11**2) * R_earth**3
        return moment


@add_common_docstring(**_frame_parameters())
class Geomagnetic(BaseMagnetic):
    """
    A coordinate or frame in the Geomagnetic (MAG) system.

    - The origin is the center of the Earth.
    - The Z-axis (+90 degrees latitude) is aligned with the Earth's magnetic north
      pole.
    - The X-axis (0 degrees longitude and 0 degrees latitude) is aligned with the
      component of the Earth's geographic north pole that is perpendicular to the
      Z-axis.

    Parameters
    ----------
    {data}
    {lonlat}
    {distance_earth}
    {magnetic_model}
    {common}

    Notes
    -----
    The position of Earth's magnetic north pole is calculated from the first three
    coefficients of the selected IGRF model per Franz & Harper (2002). The small
    offset between dipole center and Earth center is ignored.

    References
    ----------
    {igrf_reference}
    """


@add_common_docstring(**_frame_parameters())
class SolarMagnetic(BaseMagnetic):
    """
    A coordinate or frame in the Solar Magnetic (SM) system.

    - The origin is the center of the Earth.
    - The Z-axis (+90 degrees latitude) is aligned with the Earth's magnetic north
      pole.
    - The X-axis (0 degrees longitude and 0 degrees latitude) is aligned with the
      component of the Earth-Sun line that is perpendicular to the Z-axis.

    Parameters
    ----------
    {data}
    {lonlat}
    {distance_earth}
    {magnetic_model}
    {common}

    Notes
    -----
    The position of Earth's magnetic north pole is calculated from the first three
    coefficients of the selected IGRF model per Franz & Harper (2002). The small
    offset between dipole center and Earth center is ignored.

    References
    ----------
    {igrf_reference}
    """


@add_common_docstring(**_frame_parameters())
class GeocentricSolarMagnetospheric(BaseMagnetic):
    """
    A coordinate or frame in the GeocentricSolarMagnetospheric (GSM) system.

    - The origin is the center of the Earth.
    - The X-axis (0 degrees longitude and 0 degrees latitude) is aligned with the
      Earth-Sun line.
    - The Z-axis (+90 degrees latitude) is aligned with the component of the Earth's
      magnetic north pole that is perpendicular to the X-axis.

    Parameters
    ----------
    {data}
    {lonlat}
    {distance_earth}
    {magnetic_model}
    {common}

    Notes
    -----
    The position of Earth's magnetic north pole is calculated from the first three
    coefficients of the selected IGRF model per Franz & Harper (2002). The small
    offset between dipole center and Earth center is ignored.

    References
    ----------
    {igrf_reference}
    """
