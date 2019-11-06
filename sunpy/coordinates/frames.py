"""
Common solar physics coordinate systems.

This submodule implements various solar physics coordinate frames for use with
the `astropy.coordinates` module.
"""
import numpy as np

import astropy.units as u
from astropy.coordinates import Attribute, ConvertError
from astropy.coordinates.baseframe import BaseCoordinateFrame, RepresentationMapping
from astropy.coordinates.representation import (CartesianRepresentation, SphericalRepresentation,
                                                CylindricalRepresentation,
                                                UnitSphericalRepresentation)

from sunpy.sun.constants import radius as _RSUN

from .frameattributes import TimeFrameAttributeSunPy, ObserverCoordinateAttribute

from sunpy.util.decorators import add_common_docstring
from sunpy.time.time import _variables_for_parse_time_docstring

__all__ = ['HeliographicStonyhurst', 'HeliographicCarrington',
           'Heliocentric', 'Helioprojective']


def _frame_parameters():
    """
    Returns formatting dictionary to use with add_common_docstring to populate frame docstrings
    """
    ret = {}

    # Each text block is missing the first indent because it already exists in the frame docstring
    ret['data'] = ("data : `~astropy.coordinates.BaseRepresentation` or ``None``\n"
                   "        A representation object or ``None`` to have no data\n"
                   "        (or use the coordinate component arguments, see below).")
    ret['common'] = (f"obstime : {_variables_for_parse_time_docstring()['parse_time_types']}\n"
                     "        The time of the observation.  This is used to determine the\n"
                     "        position of solar-system bodies (e.g., the Sun and the Earth) as\n"
                     "        needed to define the origin and orientation of the frame.\n"
                     "    representation_type : `~astropy.coordinates.BaseRepresentation`, str, optional\n"
                     "        A representation class or string name of a representation class.\n"
                     "        This may change the valid coordinate component arguments from the\n"
                     "        defaults (see above). For example, passing\n"
                     "        ``representation_type='cartesian'`` will make the frame expect\n"
                     "        Cartesian coordinate component arguments (typically, ``x``, ``y``,\n"
                     "        and ``z``).\n"
                     "    copy : bool, optional\n"
                     "        If `True` (default), make copies of the input coordinate arrays.")
    ret['lonlat'] = ("lon : `~astropy.coordinates.Angle` or `~astropy.units.Quantity`, optional\n"
                     "        The longitude coordinate for this object (``lat`` must also be\n"
                     "        given and ``data`` must be ``None``).\n"
                     "        Not needed if ``data`` is given.\n"
                     "    lat : `~astropy.coordinates.Angle` or `~astropy.units.Quantity`, optional\n"
                     "        The latitude coordinate for this object (``lon`` must also be\n"
                     "        given and ``data`` must be ``None``).\n"
                     "        Not needed if ``data`` is given.")
    ret['radius'] = ("radius : `~astropy.units.Quantity`, optional\n"
                     "        The radial distance coordinate from Sun center for this object.\n"
                     "        Defaults to the radius of the Sun. Not needed if ``data`` is given.")
    ret['xyz'] = ("x : `~astropy.units.Quantity`, optional\n"
                  "        X-axis coordinate for this object. Not needed if ``data`` is given.\n"
                  "    y : `~astropy.units.Quantity`, optional\n"
                  "        Y-axis coordinate for this object. Not needed if ``data`` is given.\n"
                  "    z : `~astropy.units.Quantity`, optional\n"
                  "        Z-axis coordinate for this object. Not needed if ``data`` is given.")
    ret['observer'] = ("observer : `~sunpy.coordinates.frames.HeliographicStonyhurst`, str\n"
                       "        The location of the observer. If a string is provided,\n"
                       "        it must be a solar system body that can be parsed by\n"
                       "        `~sunpy.coordinates.ephemeris.get_body_heliographic_stonyhurst`\n"
                       "        at the time ``obstime``. Defaults to Earth center.")

    return ret


class SunPyBaseCoordinateFrame(BaseCoordinateFrame):
    """
    * Defines the frame attribute ``obstime`` for observation time.
    * Defines a default longitude wrap angle of 180 degrees, which can be overridden via the class
      variable ``_wrap_angle``.
    * Inject a nice way of representing the object which the coordinate represents.
    """
    obstime = TimeFrameAttributeSunPy()

    _wrap_angle = 180*u.deg

    def __init__(self, *args, **kwargs):
        self.object_name = None

        # If wrap_longitude=False is passed in, do not impose a specific wrap angle for the frame
        if not kwargs.pop('wrap_longitude', True):
            self._wrap_angle = None

        super().__init__(*args, **kwargs)

        # If obstime is specified, treat the default observer (Earth) as explicitly set
        if self.obstime is not None and self.is_frame_attr_default('observer'):
            self._attr_names_with_defaults.remove('observer')

        return

    def represent_as(self, base, s='base', in_frame_units=False):
        """
        If a frame wrap angle is set, use that wrap angle for any spherical representations.
        """
        data = super().represent_as(base, s, in_frame_units=in_frame_units)
        if self._wrap_angle is not None and \
           isinstance(data, (UnitSphericalRepresentation, SphericalRepresentation)):
            data.lon.wrap_angle = self._wrap_angle
        return data

    def __str__(self):
        """
        We override this here so that when you print a SkyCoord it shows the
        observer as the string and not the whole massive coordinate.
        """
        if getattr(self, "object_name", None):
            return f"<{self.__class__.__name__} Coordinate for '{self.object_name}'>"
        else:
            return super().__str__()


@add_common_docstring(**_frame_parameters())
class HeliographicStonyhurst(SunPyBaseCoordinateFrame):
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
    <SkyCoord (HeliographicStonyhurst: obstime=2010-01-01T00:00:45.000): (lon, lat, radius) in (deg, deg, km)
        (1., 1., 2.)>
    >>> sc.frame
    <HeliographicStonyhurst Coordinate (obstime=2010-01-01T00:00:45.000): (lon, lat, radius) in (deg, deg, km)
        (1., 1., 2.)>
    >>> sc = SkyCoord(HeliographicStonyhurst(-10*u.deg, 2*u.deg))
    >>> sc
    <SkyCoord (HeliographicStonyhurst: obstime=None): (lon, lat, radius) in (deg, deg, km)
        (-10., 2., 695700.)>
    >>> sc = SkyCoord(CartesianRepresentation(0*u.km, 45*u.km, 2*u.km),
    ...               obstime="2011/01/05T00:00:50",
    ...               frame="heliographic_stonyhurst")
    >>> sc
    <SkyCoord (HeliographicStonyhurst: obstime=2011-01-05T00:00:50.000): (lon, lat, radius) in (deg, deg, km)
    (90., 2.54480438, 45.04442252)>

    Notes
    -----
    This frame will always be converted a 3D frame where the radius defaults to
    ``rsun``.
    """
    name = "heliographic_stonyhurst"
    default_representation = SphericalRepresentation

    frame_specific_representation_info = {
        SphericalRepresentation: [RepresentationMapping(reprname='lon',
                                                        framename='lon',
                                                        defaultunit=u.deg),
                                  RepresentationMapping(reprname='lat',
                                                        framename='lat',
                                                        defaultunit=u.deg),
                                  RepresentationMapping(reprname='distance',
                                                        framename='radius',
                                                        defaultunit=None)],
        CartesianRepresentation: [RepresentationMapping(reprname='x',
                                                        framename='x'),
                                  RepresentationMapping(reprname='y',
                                                        framename='y'),
                                  RepresentationMapping(reprname='z',
                                                        framename='z')]
    }

    def __init__(self, *args, **kwargs):
        _rep_kwarg = kwargs.get('representation_type', None)

        super().__init__(*args, **kwargs)

        # Make 3D if specified as 2D
        if (self._data is not None and self._data.norm().unit is u.one
            and u.allclose(self._data.norm(), 1*u.one)):

            self._data *= _RSUN.to(u.km)


@add_common_docstring(**_frame_parameters())
class HeliographicCarrington(HeliographicStonyhurst):
    """
    A coordinate or frame in the Carrington Heliographic (HGC) system.

    - The origin is the center of the Sun.
    - The Z-axis (+90 degrees latitude) is aligned with the Sun's north pole.
    - The X-axis and Y-axis rotate with a period of 25.38 days.

    This system differs from Stonyhurst Heliographic (HGS) in its definition of longitude.

    A new instance can be created using the following signatures
    (note that if supplied, ``obstime`` must be a keyword argument)::

        HeliographicCarrington(lon, lat, obstime=obstime)
        HeliographicCarrington(lon, lat, radius, obstime=obstime)

    Parameters
    ----------
    {data}
    {lonlat}
    {radius}
    {common}

    Examples
    --------
    >>> from astropy.coordinates import SkyCoord
    >>> import sunpy.coordinates
    >>> import astropy.units as u
    >>> sc = SkyCoord(1*u.deg, 2*u.deg, 3*u.km,
    ...               frame="heliographic_carrington",
    ...               obstime="2010/01/01T00:00:30")
    >>> sc
    <SkyCoord (HeliographicCarrington: obstime=2010-01-01T00:00:30.000): (lon, lat, radius) in (deg, deg, km)
        (1., 2., 3.)>

    >>> sc = SkyCoord([1,2,3]*u.deg, [4,5,6]*u.deg, [5,6,7]*u.km,
    ...               obstime="2010/01/01T00:00:45", frame="heliographic_carrington")
    >>> sc
    <SkyCoord (HeliographicCarrington: obstime=2010-01-01T00:00:45.000): (lon, lat, radius) in (deg, deg, km)
        [(1., 4., 5.), (2., 5., 6.), (3., 6., 7.)]>

    >>> sc = SkyCoord(CartesianRepresentation(0*u.km, 45*u.km, 2*u.km),
    ...               obstime="2011/01/05T00:00:50",
    ...               frame="heliographic_carrington")
    >>> sc
    <SkyCoord (HeliographicCarrington: obstime=2011-01-05T00:00:50.000): (lon, lat, radius) in (deg, deg, km)
        (90., 2.54480438, 45.04442252)>
    """
    name = "heliographic_carrington"
    _wrap_angle = 360*u.deg


@add_common_docstring(**_frame_parameters())
class Heliocentric(SunPyBaseCoordinateFrame):
    """
    A coordinate or frame in the Heliocentric system, which is observer-based.

    - The origin is the center of the Sun.
    - The Z-axis is aligned with the Sun-observer line.
    - The Y-axis is aligned with the component of the vector to the Sun's north pole that is
      perpendicular to the Z-axis.

    This frame defaults to a Cartesian component representation, which is known as Heliocentric
    Cartesian (HCC).  This frame can also be represented using cylindrical components, where
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
    ...               obstime="2011/01/05T00:00:50", frame="heliocentric")
    >>> sc
    <SkyCoord (Heliocentric: obstime=2011-01-05T00:00:50.000, observer=<HeliographicStonyhurst Coordinate for 'earth'>): (x, y, z) in km
        (10., 1., 2.)>

    >>> sc = SkyCoord([1,2]*u.km, [3,4]*u.m, [5,6]*u.cm, frame="heliocentric", obstime="2011/01/01T00:00:54")

    >>> sc
    <SkyCoord (Heliocentric: obstime=2011-01-01T00:00:54.000, observer=<HeliographicStonyhurst Coordinate for 'earth'>): (x, y, z) in (km, m, cm)
        [(1., 3., 5.), (2., 4., 6.)]>

    >>> sc = SkyCoord(CylindricalRepresentation(10*u.km, 60*u.deg, 10*u.km),
    ...               obstime="2011/01/05T00:00:50",
    ...               frame="heliocentric")
    >>> sc
    <SkyCoord (Heliocentric: obstime=2011-01-05T00:00:50.000, observer=<HeliographicStonyhurst Coordinate for 'earth'>): (x, y, z) in km
        (5., 8.66025404, 10.)>
    """
    default_representation = CartesianRepresentation

    frame_specific_representation_info = {
        CylindricalRepresentation: [RepresentationMapping('phi', 'psi', u.deg)]
    }

    observer = ObserverCoordinateAttribute(HeliographicStonyhurst, default="earth")


@add_common_docstring(**_frame_parameters())
class Helioprojective(SunPyBaseCoordinateFrame):
    """
    A coordinate or frame in the Helioprojective Cartesian (HPC) system, which is observer-based.

    - The origin is the location of the observer.
    - ``theta_x`` is the angle relative to the plane containing the Sun-observer line and the Sun's
      rotation axis, with positive values in the direction of the Sun's west limb.
    - ``theta_y`` is the angle relative to the Sun's equatorial plane, with positive values in the
      direction of the Sun's north pole.
    - ``distance`` is the Sun-observer distance.

    This system is frequently used in a projective form without ``distance`` specified.  For
    observations looking very close to the center of the Sun, where the small-angle approximation
    is appropriate, ``theta_x`` and ``theta_y`` can be approximated as Cartesian components.

    A new instance can be created using the following signatures
    (note that if supplied, ``obstime`` and ``observer`` must be keyword arguments)::

        Helioprojective(theta_x, theta_y, obstime=obstime, observer=observer)
        Helioprojective(theta_x, theta_y, distance, obstime=obstime, observer=observer)

    Parameters
    ----------
    {data}
    Tx : `~astropy.coordinates.Angle` or `~astropy.units.Quantity`
        The theta_x coordinate for this object. Not needed if ``data`` is given.
    Ty : `~astropy.coordinates.Angle` or `~astropy.units.Quantity`
        The theta_y coordinate for this object. Not needed if ``data`` is given.
    distance : `~astropy.units.Quantity`
        The distance coordinate from the observer for this object.
        Not needed if ``data`` is given.
    {observer}
    rsun : `~astropy.units.Quantity`
        The physical (length) radius of the Sun. Used to calculate the position
        of the limb for calculating distance from the observer to the
        coordinate. Defaults to the solar radius.
    {common}

    Examples
    --------
    >>> from astropy.coordinates import SkyCoord
    >>> import sunpy.coordinates
    >>> import astropy.units as u
    >>> sc = SkyCoord(0*u.deg, 0*u.deg, 5*u.km, obstime="2010/01/01T00:00:00",
    ...               frame="helioprojective")
    >>> sc
    <SkyCoord (Helioprojective: obstime=2010-01-01T00:00:00.000, rsun=695700.0 km, observer=<HeliographicStonyhurst Coordinate for 'earth'>): (Tx, Ty, distance) in (arcsec, arcsec, km)
        (0., 0., 5.)>
    >>> sc = SkyCoord(0*u.deg, 0*u.deg, obstime="2010/01/01T00:00:00", frame="helioprojective")
    >>> sc
    <SkyCoord (Helioprojective: obstime=2010-01-01T00:00:00.000, rsun=695700.0 km, observer=<HeliographicStonyhurst Coordinate for 'earth'>): (Tx, Ty) in arcsec
        (0., 0.)>
    >>> sc = SkyCoord(CartesianRepresentation(1*u.AU, 1e5*u.km, -2e5*u.km),
    ...               obstime="2011/01/05T00:00:50",
    ...               frame="helioprojective")
    >>> sc
    <SkyCoord (Helioprojective: obstime=2011-01-05T00:00:50.000, rsun=695700.0 km, observer=<HeliographicStonyhurst Coordinate for 'earth'>): (Tx, Ty, distance) in (arcsec, arcsec, AU)
        (137.87948623, -275.75878762, 1.00000112)>
    """
    default_representation = SphericalRepresentation

    frame_specific_representation_info = {
        SphericalRepresentation: [RepresentationMapping(reprname='lon',
                                                        framename='Tx',
                                                        defaultunit=u.arcsec),
                                  RepresentationMapping(reprname='lat',
                                                        framename='Ty',
                                                        defaultunit=u.arcsec),
                                  RepresentationMapping(reprname='distance',
                                                        framename='distance',
                                                        defaultunit=None)],

        UnitSphericalRepresentation: [RepresentationMapping(reprname='lon',
                                                            framename='Tx',
                                                            defaultunit=u.arcsec),
                                      RepresentationMapping(reprname='lat',
                                                            framename='Ty',
                                                            defaultunit=u.arcsec)],
    }

    rsun = Attribute(default=_RSUN.to(u.km))
    observer = ObserverCoordinateAttribute(HeliographicStonyhurst, default="earth")

    def calculate_distance(self):
        """
        This method calculates the third coordinate of the Helioprojective
        frame. It assumes that the coordinate point is on the disk of the Sun
        at the rsun radius.

        If a point in the frame is off limb then NaN will be returned.

        Returns
        -------
        new_frame : `~sunpy.coordinates.frames.HelioProjective`
            A new frame instance with all the attributes of the original but
            now with a third coordinate.
        """
        # Skip if we already are 3D
        distance = self.spherical.distance
        if not (distance.unit is u.one and u.allclose(distance, 1*u.one)):
            return self

        if not isinstance(self.observer, BaseCoordinateFrame):
            raise ConvertError("Cannot calculate distance to the solar disk "
                               "for observer '{}' "
                               "without `obstime` being specified.".format(self.observer))

        rep = self.represent_as(UnitSphericalRepresentation)
        lat, lon = rep.lat, rep.lon
        alpha = np.arccos(np.cos(lat) * np.cos(lon)).to(lat.unit)
        c = self.observer.radius**2 - self.rsun**2
        b = -2 * self.observer.radius * np.cos(alpha)
        # Ingore sqrt of NaNs
        with np.errstate(invalid='ignore'):
            d = ((-1*b) - np.sqrt(b**2 - 4*c)) / 2

        return self.realize_frame(SphericalRepresentation(lon=lon,
                                                          lat=lat,
                                                          distance=d))
