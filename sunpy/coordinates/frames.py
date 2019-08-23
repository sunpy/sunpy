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


@add_common_docstring(**_variables_for_parse_time_docstring())
class HeliographicStonyhurst(SunPyBaseCoordinateFrame):
    """
    A coordinate or frame in the Stonyhurst Heliographic system.

    In a Cartesian representation this is also known as the Heliocentric
    Earth Equatorial (HEEQ) system. This frame has its origin at the solar
    center and the north pole above the solar north pole, and the zero line on
    longitude pointing towards the Earth. If the ``data`` parameter
    is given, the positional parameters for the coordinate frame
    (``lon``, ``lat``, ``radius``) do not need to be given.

    A new instance can be created using the following signatures
    (note that ``obstime`` and ``representation_type`` must be supplied as
    keywords)::

        HeliographicStonyhurst(lon, lat, obstime)
        HeliographicStonyhurst(lon, lat, radius, obstime)
        HeliographicStonyhurst(x, y, z, obstime, representation_type='cartesian')

    Parameters
    ----------
    data : `~astropy.coordinates.BaseRepresentation` or `None`
        A representation object or None to have no data.
    lon : `~astropy.coordinates.Angle` or `~astropy.units.Quantity`, optional
        The longitude for this object (``lat`` must also be given and
        ``data`` must be None). Not needed if ``data`` is given.
    lat : `~astropy.coordinates.Angle` or `~astropy.units.Quantity`, optional
        The latitude for this object (``lon`` must also be given and
        ``data`` must be None). Not needed if ``data`` is given.
    radius : `~astropy.units.Quantity`, optional
        The radial distance for this object. Defaults to the solar
        radius. Not needed if ``data`` is given.
    obstime : {parse_time_types}
        The date and time of the observation.

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
    rsun.
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
        # If representation was explicitly passed, do not change the rep.
        if not _rep_kwarg:
            if isinstance(self._data, UnitSphericalRepresentation):
                self._data = self.spherical

    def represent_as(self, base, s='base', in_frame_units=False):
        """
        Unless the requested representation is UnitSphericalRepresentation, scale a coordinate with
        dimensionless length so that it has the length of the solar radius.
        """
        data = super().represent_as(base, s, in_frame_units=in_frame_units)
        if not isinstance(data, UnitSphericalRepresentation) and \
           data.norm().unit is u.one and u.allclose(data.norm(), 1*u.one):
            data *= _RSUN.to(u.km)
        return data


@add_common_docstring(**_variables_for_parse_time_docstring())
class HeliographicCarrington(HeliographicStonyhurst):
    """
    A coordinate or frame in the Carrington Heliographic system.

    - The origin is the centre of the Sun
    - The z-axis is aligned with the Sun's north pole
    - The x and y axes rotate with a period of 25.38 days.

    This frame differs from the Stonyhurst version in the definition of the
    longitude, which is defined using the time-dependent offset described
    above. If the ``data`` parameter is given, the positional
    parameters for the coordinate frame (``lon``, ``lat``, ``radius``) do not need to be
    given.

    Parameters
    ----------
    data : `~astropy.coordinates.BaseRepresentation` or None
        A representation object or None to have no data.
    lon : `~astropy.coordinates.Angle`, optional
        The longitude for this object (``lat`` must also be given and
        ``data`` must be None). Not needed if ``data`` is given.
    lat : `~astropy.coordinates.Angle`, optional
        The latitude for this object (``lon`` must also be given and
        ``data`` must be None). Not needed if ``data`` is given.
    radius : `~astropy.units.Quantity`, optional
        The radial distance for this object. Defaults to the solar radius.
        Not needed if ``data`` is given.
    obstime : {parse_time_types}
        The date and time of the observation.

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


@add_common_docstring(**_variables_for_parse_time_docstring())
class Heliocentric(SunPyBaseCoordinateFrame):
    """
    A coordinate or frame in the Heliocentric system.

    - The origin is the centre of the Sun
    - The z-axis points from the centre of the Sun to the observer.
    - The y-axis is perpendicular to the z-axis, and lies in the plane that
      contains the z-axis and the solar rotation axis, pointing towards the
      Sun's north pole.

    This frame may either be specified in Cartesian or cylindrical
    representation. Cylindrical representation replaces (``x``, ``y``) with (``rho``, ``psi``)
    where ``rho`` is the impact parameter and ``psi`` is the position angle.
    ``psi`` is measured relative to the west limb, rather than solar north, so is shifted
    by 90 degrees compared to the convention of the Heliocentric Radial system.
    If the ``data`` parameter is given, the positional parameters
    for the coordinate frame (``x``, ``y``, ``z``) do not need to be given.

    A new instance can be created using the following signatures
    (note that ``obstime`` and ``representation_type`` must be supplied as
    keywords)::

        Heliocentric(x, y, z, obstime)
        Heliocentric(rho, psi, z, obstime, representation_type='cylindrical')

    Parameters
    ----------
    data : `~astropy.coordinates.BaseRepresentation` or None
        A representation object. If specified, other parameters must
        be in keyword form and if x, y and z are specified, it must
        be None.
    x : `~astropy.units.Quantity`, optional
        X-axis coordinate. Not needed if ``data`` is given.
    y : `~astropy.units.Quantity`, optional
        Y-axis coordinate. Not needed if ''data'' is given.
    z : `~astropy.units.Quantity`, optional
        Z-axis coordinate. Not needed if ``data`` is given.
    observer : `~sunpy.coordinates.frames.HeliographicStonyhurst`, str
        The coordinate of the observer in the solar system. If you
        supply a string, it must be a solar system body that can be
        parsed by `~sunpy.coordinates.ephemeris.get_body_heliographic_stonyhurst`.
        Defaults to the Earth.
    obstime : {parse_time_types}
        The date and time of the observation.

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


@add_common_docstring(**_variables_for_parse_time_docstring())
class Helioprojective(SunPyBaseCoordinateFrame):
    """
    A coordinate or frame in the Helioprojective (Cartesian) system.

    This is a projective coordinate system centered around the observer.
    It is a full spherical coordinate system with position given as longitude
    theta_x and latitude theta_y. If the ``data`` parameter is given,
    the positional parameters for the coordinate frame (``Tx``, ``Ty``, ``distance``)
    do not need to be given.

    Parameters
    ----------
    data : `~astropy.coordinates.BaseRepresentation` or None
        A representation object. If specified, other parameters must
        be in keyword form.
    Tx : `~astropy.coordinates.Angle` or `~astropy.units.Quantity`
        Theta_x coordinate. Not needed if ``data`` is given.
    Ty : `~astropy.coordinates.Angle` or `~astropy.units.Quantity`
        Theta_y coordinate. Not needed if ``data`` is given.
    distance : `~astropy.units.Quantity`
        The radial distance from the observer to the coordinate point.
        Not needed if ``data`` is given.
    obstime : {parse_time_types}
        The date and time of the observation.
    observer : `~sunpy.coordinates.frames.HeliographicStonyhurst`, str
        The coordinate of the observer in the solar system. If you supply a string,
        it must be a solar system body that can be parsed by
        `~sunpy.coordinates.ephemeris.get_body_heliographic_stonyhurst`. Defaults
        to the Earth.
    rsun : `~astropy.units.Quantity`
        The physical (length) radius of the Sun. Used to calculate the position
        of the limb for calculating distance from the observer to the
        coordinate. Defaults to the solar radius.

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
