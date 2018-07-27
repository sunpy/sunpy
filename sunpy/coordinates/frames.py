"""
Common solar physics coordinate systems.

This submodule implements various solar physics coordinate frames for use with
the `astropy.coordinates` module.
"""
from __future__ import absolute_import, division

import numpy as np
from astropy import units as u
from astropy.coordinates import Attribute, ConvertError
from astropy.coordinates.baseframe import (BaseCoordinateFrame,
                                           RepresentationMapping)
from astropy.coordinates.representation import (CartesianRepresentation,
                                                CylindricalRepresentation,
                                                SphericalRepresentation,
                                                UnitSphericalRepresentation)
from sunpy import sun

from .frameattributes import (ObserverCoordinateAttribute,
                              TimeFrameAttributeSunPy)

try:
    from astropy.units import allclose as quantity_allclose
except ImportError:
    from astropy.tests.helper import quantity_allclose


RSUN_METERS = sun.constants.get('radius').si.to(u.m)
DSUN_METERS = sun.constants.get('mean distance').si.to(u.m)

__all__ = ['HeliographicStonyhurst', 'HeliographicCarrington', 'Heliocentric',
           'Helioprojective']


class HeliographicStonyhurst(BaseCoordinateFrame):
    """
    A coordinate or frame in the Stonyhurst Heliographic system.

    In a cartesian representation this is also known as the Heliocentric
    Earth Equatorial (HEEQ) system. This frame has its origin at the solar
    centre and the north pole above the solar north pole, and the zero line on
    longitude pointing towards the Earth.

    A new instance can be created using the following signatures
    (note that all the arguments must be supplied as keywords)::

        HeliographicStonyhurst(lon, lat, obstime)
        HeliographicStonyhurst(lon, lat, radius, obstime)
        HeliographicStonyhurst(x, y, z, obstime, representation='cartesian')

    Parameters
    ----------
    representation : `~astropy.coordinates.BaseRepresentation` or `None`
        A representation object or None to have no data.
    lon : `~astropy.coordinates.Angle`, optional
        The longitude for this object (``lat`` must also be given and
        ``representation`` must be None).
    lat : `~astropy.coordinates.Angle`, optional
        The latitude for this object (``lon`` must also be given and
        ``representation`` must be None).
    radius : `~astropy.units.Quantity`, optional
        This quantity holds the radial distance. If not specified, it is, by
        default, the radius of the photosphere.
    x : `~astropy.units.Quantity`, optional
        x coordinate.
    y : `~astropy.units.Quantity`, optional
        y coordinate.
    z : `~astropy.units.Quantity`, optional
        z coordinate.
    obstime: `~sunpy.time.Time`
        The date and time of the observation, used to convert to heliographic
        carrington coordinates.

    Examples
    --------
    >>> from astropy.coordinates import SkyCoord
    >>> import sunpy.coordinates
    >>> import astropy.units as u
    >>> sc = SkyCoord(1*u.deg, 1*u.deg, 2*u.km,
    ...               frame="heliographic_stonyhurst",
    ...               obstime="2010/01/01T00:00:45")
    >>> sc
    <SkyCoord (HeliographicStonyhurst: obstime=2010-01-01 00:00:45): (lon, lat, radius) in (deg, deg, km)
        (1., 1., 2.)>
    >>> sc.frame
    <HeliographicStonyhurst Coordinate (obstime=2010-01-01 00:00:45): (lon, lat, radius) in (deg, deg, km)
        (1., 1., 2.)>
    >>> sc = SkyCoord(HeliographicStonyhurst(-10*u.deg, 2*u.deg))
    >>> sc
    <SkyCoord (HeliographicStonyhurst: obstime=None): (lon, lat, radius) in (deg, deg, km)
        (-10., 2., 695508.)>

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

    obstime = TimeFrameAttributeSunPy()

    _default_wrap_angle = 180*u.deg

    def __init__(self, *args, **kwargs):
        _rep_kwarg = kwargs.get('representation', None)
        wrap = kwargs.pop('wrap_longitude', True)

        if ('radius' in kwargs and kwargs['radius'].unit is u.one and
                quantity_allclose(kwargs['radius'], 1*u.one)):
            kwargs['radius'] = RSUN_METERS.to(u.km)

        super(HeliographicStonyhurst, self).__init__(*args, **kwargs)

        # Make 3D if specified as 2D
        # If representation was explicitly passed, do not change the rep.
        if not _rep_kwarg:
            # If we were passed a 3D rep extract the distance, otherwise
            # calculate it from RSUN.
            if isinstance(self._data, UnitSphericalRepresentation):
                distance = RSUN_METERS.to(u.km)
                self._data = SphericalRepresentation(lat=self._data.lat,
                                                     lon=self._data.lon,
                                                     distance=distance)

        if wrap and isinstance(self._data, (UnitSphericalRepresentation, SphericalRepresentation)):
            self._data.lon.wrap_angle = self._default_wrap_angle


class HeliographicCarrington(HeliographicStonyhurst):
    """
    A coordinate or frame in the Carrington Heliographic system.

    - The origin is the centre of the Sun
    - The z-axis is aligned with the Sun's north pole
    - The x and y axes rotate with a period of 25.38 days. The line of zero
      longitude passed through the disk centre as seen from Earth at
      21:36 on 9th Nov 1853.

    This frame differs from the Stonyhurst version in the definition of the
    longitude, which is defined using the time-dependant offset described
    above.

    Parameters
    ----------
    representation: `~astropy.coordinates.BaseRepresentation` or None.
        A representation object. If specified, other parameters must
        be in keyword form.
    lon: `Angle` object.
        The longitude for this object (``lat`` must also be given and
        ``representation`` must be None).
    lat: `Angle` object.
        The latitude for this object (``lon`` must also be given and
        ``representation`` must be None).
    radius: `astropy.units.Quantity` object, optional, must be keyword.
        This quantity holds the radial distance. Defaults to the solar radius.
    obstime: SunPy Time
        The date and time of the observation, used to convert to heliographic
        carrington coordinates.

    Examples
    --------
    >>> from astropy.coordinates import SkyCoord
    >>> import sunpy.coordinates
    >>> import astropy.units as u
    >>> sc = SkyCoord(1*u.deg, 2*u.deg, 3*u.km,
    ...               frame="heliographic_carrington",
    ...               obstime="2010/01/01T00:00:30")
    >>> sc
    <SkyCoord (HeliographicCarrington: obstime=2010-01-01 00:00:30): (lon, lat, radius) in (deg, deg, km)
        (1., 2., 3.)>

    >>> sc = SkyCoord([1,2,3]*u.deg, [4,5,6]*u.deg, [5,6,7]*u.km,
    ...               obstime="2010/01/01T00:00:45", frame="heliographic_carrington")
    >>> sc
    <SkyCoord (HeliographicCarrington: obstime=2010-01-01 00:00:45): (lon, lat, radius) in (deg, deg, km)
        [(1., 4., 5.), (2., 5., 6.), (3., 6., 7.)]>
    """

    name = "heliographic_carrington"
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

        UnitSphericalRepresentation: [RepresentationMapping(reprname='lon',
                                                            framename='lon',
                                                            defaultunit=u.deg),
                                      RepresentationMapping(reprname='lat',
                                                            framename='lat',
                                                            defaultunit=u.deg)],
    }

    _default_wrap_angle = 360*u.deg
    obstime = TimeFrameAttributeSunPy()


class Heliocentric(BaseCoordinateFrame):
    """
    A coordinate or frame in the Heliocentric system.

    - The origin is the centre of the Sun
    - The z-axis points from the centre of the Sun to the observer.
    - The y-axis is perpendicular to the z-axis, and lies in the plane that
      contains the z-axis and the solar rotation axis, pointing towards the
      Sun's north pole.

    This frame may either be specified in Cartesian or cylindrical
    representation. Cylindrical representation replaces (x, y) with (rho, psi)
    where rho is the impact parameter and psi is the position angle in degrees.

    Parameters
    ----------
    representation: `~astropy.coordinates.BaseRepresentation` or None.
        A representation object. If specified, other parameters must
        be in keyword form and if x, y and z are specified, it must
        be None.
    x: `Quantity` object.
        X-axis coordinate, optional, must be keyword.
    y: `Quantity` object.
        Y-axis coordinate, optional, must be keyword.
    z: `Quantity` object. Shared by both representations.
        Z-axis coordinate, optional, must be keyword.
    observer: `~sunpy.coordinates.frames.HeliographicStonyhurst`, optional
        The coordinate of the observer in the solar system. Defaults to the
        Earth.
    obstime: SunPy Time
        The date and time of the observation, used to convert to heliographic
        carrington coordinates.

    Examples
    --------

    >>> from astropy.coordinates import SkyCoord, CartesianRepresentation
    >>> import sunpy.coordinates
    >>> import astropy.units as u

    >>> sc = SkyCoord(CartesianRepresentation(10*u.km, 1*u.km, 2*u.km),
    ...               obstime="2011/01/05T00:00:50", frame="heliocentric")
    >>> sc
    <SkyCoord (Heliocentric: obstime=2011-01-05 00:00:50, observer=<HeliographicStonyhurst Coordinate (obstime=2011-01-05 00:00:50): (lon, lat, radius) in (deg, deg, AU)
        (0., -3.43939103, 0.98334411)>): (x, y, z) in km
        (10., 1., 2.)>

    >>> sc = SkyCoord([1,2]*u.km, [3,4]*u.m, [5,6]*u.cm, frame="heliocentric", obstime="2011/01/01T00:00:54")

    >>> sc
    <SkyCoord (Heliocentric: obstime=2011-01-01 00:00:54, observer=<HeliographicStonyhurst Coordinate (obstime=2011-01-01 00:00:54): (lon, lat, radius) in (deg, deg, AU)
        (0., -2.97725356, 0.98335586)>): (x, y, z) in (km, m, cm)
        [(1., 3., 5.), (2., 4., 6.)]>
    """

    default_representation = CartesianRepresentation

    _frame_specific_representation_info = {
        CylindricalRepresentation: [RepresentationMapping('phi', 'psi', u.deg)]
    }

    obstime = TimeFrameAttributeSunPy()
    observer = ObserverCoordinateAttribute(HeliographicStonyhurst, default="earth")


class Helioprojective(BaseCoordinateFrame):
    """
    A coordinate or frame in the Helioprojective (Cartesian) system.

    This is a projective coordinate system centered around the observer.
    It is a full spherical coordinate system with position given as longitude
    theta_x and latitude theta_y.

    Parameters
    ----------
    representation: `~astropy.coordinates.BaseRepresentation` or None.
        A representation object. If specified, other parameters must
        be in keyword form.
    Tx: `~astropy.coordinates.Angle`  or `~astropy.units.Quantity`
        X-axis coordinate.
    Ty: `~astropy.coordinates.Angle`  or `~astropy.units.Quantity`
        Y-axis coordinate.
    distance: `~astropy.units.Quantity`
        The radial distance from the observer to the coordinate point.
    obstime: SunPy Time
        The date and time of the observation, used to convert to heliographic
        carrington coordinates.
    observer: `~sunpy.coordinates.frames.HeliographicStonyhurst`
        The coordinate of the observer in the solar system.
    rsun: `~astropy.units.Quantity`
        The physical (length) radius of the Sun. Used to calculate the position
        of the limb for calculating distance from the observer to the
        coordinate.

    Examples
    --------
    >>> from astropy.coordinates import SkyCoord
    >>> import sunpy.coordinates
    >>> import astropy.units as u
    >>> sc = SkyCoord(0*u.deg, 0*u.deg, 5*u.km, obstime="2010/01/01T00:00:00",
    ...               frame="helioprojective")
    >>> sc
    <SkyCoord (Helioprojective: obstime=2010-01-01 00:00:00, rsun=695508.0 km, observer=<HeliographicStonyhurst Coordinate (obstime=2010-01-01 00:00:00): (lon, lat, radius) in (deg, deg, AU)
        (0., -3.00724817, 0.98330294)>): (Tx, Ty, distance) in (arcsec, arcsec, km)
        (0., 0., 5.)>
    >>> sc = SkyCoord(0*u.deg, 0*u.deg, obstime="2010/01/01T00:00:00", frame="helioprojective")
    >>> sc
    <SkyCoord (Helioprojective: obstime=2010-01-01 00:00:00, rsun=695508.0 km, observer=<HeliographicStonyhurst Coordinate (obstime=2010-01-01 00:00:00): (lon, lat, radius) in (deg, deg, AU)
        (0., -3.00724817, 0.98330294)>): (Tx, Ty) in arcsec
        (0., 0.)>
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

    obstime = TimeFrameAttributeSunPy()
    rsun = Attribute(default=RSUN_METERS.to(u.km))
    observer = ObserverCoordinateAttribute(HeliographicStonyhurst, default="earth")

    def __init__(self, *args, **kwargs):
        wrap = kwargs.pop('wrap_longitude', True)

        BaseCoordinateFrame.__init__(self, *args, **kwargs)

        if wrap and isinstance(self._data, (UnitSphericalRepresentation, SphericalRepresentation)):
            self._data.lon.wrap_angle = 180*u.deg

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
        if (isinstance(self._data, SphericalRepresentation) and
                not (self.distance.unit is u.one and quantity_allclose(self.distance, 1*u.one))):
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
        d = ((-1*b) - np.sqrt(b**2 - 4*c)) / 2

        return self.realize_frame(SphericalRepresentation(lon=lon,
                                                          lat=lat,
                                                          distance=d))
