"""
Common solar physics coordinate systems.

This submodule implements various solar physics coordinate frames for use with
the `astropy.coordinates` module.
"""
from __future__ import absolute_import, division

import numpy as np

from astropy import units as u
from astropy.coordinates.representation import (CartesianRepresentation,
                                                UnitSphericalRepresentation,
                                                SphericalRepresentation)
from astropy.coordinates.baseframe import (BaseCoordinateFrame,
                                           RepresentationMapping)
from astropy.coordinates import Attribute, ConvertError

from sunpy import sun
from .representation import (SphericalWrap180Representation, UnitSphericalWrap180Representation)

from .represenations import SouthPoleSphericalRepresentation, UnitSouthPoleSphericalRepresentation
from .frameattributes import TimeFrameAttributeSunPy, ObserverCoordinateAttribute

RSUN_METERS = sun.constants.get('radius').si.to(u.m)
DSUN_METERS = sun.constants.get('mean distance').si.to(u.m)

__all__ = ['HeliographicStonyhurst', 'HeliographicCarrington', 'Heliocentric',
           'Helioprojective', 'HelioprojectiveRadial']


class HeliographicStonyhurst(BaseCoordinateFrame):
    """
    A coordinate or frame in the Stonyhurst Heliographic
    system.

    This frame has its origin at the solar centre and the north pole above the
    solar north pole, and the zero line on longitude pointing towards the
    Earth.

    Parameters
    ----------
    representation: `~astropy.coordinates.BaseRepresentation` or `None`
        A representation object or None to have no data.
    lon: `Angle` object.
        The longitude for this object (``lat`` must also be given and
        ``representation`` must be None).
    lat: `Angle` object.
        The latitude for this object (``lon`` must also be given and
        ``representation`` must be None).
    radius: `astropy.units.Quantity` object.
        This quantity holds the radial distance. If not specified, it is, by
        default, the radius of the photosphere. Optional.
    obstime: SunPy Time
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
    >>> sc # doctest: +FLOAT_CMP
    <SkyCoord (HeliographicStonyhurst: obstime=2010-01-01 00:00:45): (lon, lat, radius) in (deg, deg, km)
        ( 1.,  1.,  2.)>
    >>> sc.frame # doctest: +FLOAT_CMP
    <HeliographicStonyhurst Coordinate (obstime=2010-01-01 00:00:45): (lon, lat, radius) in (deg, deg, km)
        ( 1.,  1.,  2.)>
    >>> sc = SkyCoord(HeliographicStonyhurst(-10*u.deg, 2*u.deg))
    >>> sc # doctest: +FLOAT_CMP
    <SkyCoord (HeliographicStonyhurst: obstime=None): (lon, lat, radius) in (deg, deg, km)
        (-10.,  2.,  695508.)>

    Notes
    -----
    This frame will always be converted a 3D frame where the radius defaults to
    rsun.
    """
    name = "heliographic_stonyhurst"
    default_representation = SphericalWrap180Representation

    _frame_specific_representation_info = {
        'spherical': [
            RepresentationMapping('lon', 'lon', 'recommended'),
            RepresentationMapping('lat', 'lat', 'recommended'),
            RepresentationMapping('distance', 'radius', 'recommended')
        ],
        'sphericalwrap180': [
            RepresentationMapping('lon', 'lon', 'recommended'),
            RepresentationMapping('lat', 'lat', 'recommended'),
            RepresentationMapping('distance', 'radius', 'recommended')
        ]
    }

    obstime = TimeFrameAttributeSunPy()

    def __init__(self, *args, **kwargs):
        _rep_kwarg = kwargs.get('representation', None)

        super(HeliographicStonyhurst, self).__init__(*args, **kwargs)

        # Make 3D if specified as 2D
        # If representation was explicitly passed, do not change the rep.
        if not _rep_kwarg:
            # If we were passed a 3D rep extract the distance, otherwise
            # calculate it from RSUN.
            distance = None
            if isinstance(self._data, SphericalRepresentation):
                distance = self._data.distance
            elif isinstance(self._data, UnitSphericalRepresentation):
                distance = RSUN_METERS.to(u.km)

            if distance is not None:
                self._data = self.default_representation(lat=self._data.lat,
                                                         lon=self._data.lon,
                                                         distance=distance)


class HeliographicCarrington(HeliographicStonyhurst):
    """
    A coordinate or frame in the Carrington Heliographic
    system.
    This frame differs from the Stonyhurst version in the
    definition of the longitude, which is defined using
    an offset which is a time-dependent scalar value.

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
        This quantity holds the radial distance. If not specified, it is, by
        default, the solar radius. Optional, must be keyword.
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
    >>> sc # doctest: +FLOAT_CMP
    <SkyCoord (HeliographicCarrington: obstime=2010-01-01 00:00:30): (lon, lat, radius) in (deg, deg, km)
        ( 1.,  2.,  3.)>

    >>> sc = SkyCoord([1,2,3]*u.deg, [4,5,6]*u.deg, [5,6,7]*u.km,
    ...               obstime="2010/01/01T00:00:45", frame="heliographic_carrington")
    >>> sc # doctest: +FLOAT_CMP
    <SkyCoord (HeliographicCarrington: obstime=2010-01-01 00:00:45): (lon, lat, radius) in (deg, deg, km)
        [( 1.,  4.,  5.), ( 2.,  5.,  6.), ( 3.,  6.,  7.)]>
    """

    name = "heliographic_carrington"
    default_representation = SphericalRepresentation

    _frame_specific_representation_info = {
        'spherical': [
            RepresentationMapping('lon', 'lon', 'recommended'),
            RepresentationMapping('lat', 'lat', 'recommended'),
            RepresentationMapping('distance', 'radius', 'recommended')
        ],
        'sphericalwrap180': [
            RepresentationMapping('lon', 'lon', 'recommended'),
            RepresentationMapping('lat', 'lat', 'recommended'),
            RepresentationMapping('distance', 'radius', 'recommended')
        ]
    }

    obstime = TimeFrameAttributeSunPy()


class Heliocentric(BaseCoordinateFrame):
    """
    A coordinate or frame in the Heliocentric system.
    This frame may either be specified in Cartesian
    or cylindrical representation.
    Cylindrical representation replaces (x, y) with
    (rho, psi) where rho is the impact parameter and
    psi is the position angle in degrees.

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
    observer: `~sunpy.coordinates.frames.HeliographicStonyhurst`
        The coordinate of the observer in the solar system.
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
    >>> sc # doctest: +FLOAT_CMP
    <SkyCoord (Heliocentric: obstime=2011-01-05 00:00:50, observer=<HeliographicStonyhurst Coordinate (obstime=2011-01-05 00:00:50): (lon, lat, radius) in (deg, deg, AU)
        ( 0., -3.43939103,  0.98334411)>): (x, y, z) in km
        ( 10.,  1.,  2.)>

    >>> sc = SkyCoord([1,2]*u.km, [3,4]*u.m, [5,6]*u.cm, frame="heliocentric", obstime="2011/01/01T00:00:54")

    >>> sc # doctest: +FLOAT_CMP
    <SkyCoord (Heliocentric: obstime=2011-01-01 00:00:54, observer=<HeliographicStonyhurst Coordinate (obstime=2011-01-01 00:00:54): (lon, lat, radius) in (deg, deg, AU)
        ( 0., -2.97725356,  0.98335586)>): (x, y, z) in (km, m, cm)
        [( 1.,  3.,  5.), ( 2.,  4.,  6.)]>
    """

    default_representation = CartesianRepresentation

    _frame_specific_representation_info = {
        'cylindrical': [RepresentationMapping('phi', 'psi', u.deg)]
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
    >>> sc # doctest: +FLOAT_CMP
    <SkyCoord (Helioprojective: obstime=2010-01-01 00:00:00, rsun=695508.0 km, observer=<HeliographicStonyhurst Coordinate (obstime=2010-01-01 00:00:00): (lon, lat, radius) in (deg, deg, AU)
        ( 0., -3.00724817,  0.98330294)>): (Tx, Ty, distance) in (arcsec, arcsec, km)
        ( 0.,  0.,  5.)>
    >>> sc = SkyCoord(0*u.deg, 0*u.deg, obstime="2010/01/01T00:00:00", frame="helioprojective")
    >>> sc # doctest: +FLOAT_CMP
    <SkyCoord (Helioprojective: obstime=2010-01-01 00:00:00, rsun=695508.0 km, observer=<HeliographicStonyhurst Coordinate (obstime=2010-01-01 00:00:00): (lon, lat, radius) in (deg, deg, AU)
        ( 0., -3.00724817,  0.98330294)>): (Tx, Ty) in arcsec
        ( 0.,  0.)>
    """

    default_representation = SphericalWrap180Representation

    _frame_specific_representation_info = {
        'spherical': [
            RepresentationMapping('lon', 'Tx', u.arcsec),
            RepresentationMapping('lat', 'Ty', u.arcsec),
            RepresentationMapping('distance', 'distance', u.km)
        ],
        'sphericalwrap180': [
            RepresentationMapping('lon', 'Tx', u.arcsec),
            RepresentationMapping('lat', 'Ty', u.arcsec),
            RepresentationMapping('distance', 'distance', u.km)
        ],
        'unitspherical': [
            RepresentationMapping('lon', 'Tx', u.arcsec),
            RepresentationMapping('lat', 'Ty', u.arcsec)
        ],
        'unitsphericalwrap180': [
            RepresentationMapping('lon', 'Tx', u.arcsec),
            RepresentationMapping('lat', 'Ty', u.arcsec)
        ]
    }

    obstime = TimeFrameAttributeSunPy()
    rsun = Attribute(default=RSUN_METERS.to(u.km))
    observer = ObserverCoordinateAttribute(HeliographicStonyhurst, default="earth")

    def __init__(self, *args, **kwargs):
        _rep_kwarg = kwargs.get('representation', None)

        BaseCoordinateFrame.__init__(self, *args, **kwargs)

        # Convert from Spherical to SphericalWrap180
        # If representation was explicitly passed, do not change the rep.
        if not _rep_kwarg:
            # The base __init__ will make this a UnitSphericalRepresentation
            # This makes it Wrap180 instead
            if isinstance(self._data, UnitSphericalRepresentation):
                self._data = UnitSphericalWrap180Representation(
                    lat=self._data.lat, lon=self._data.lon)
                self.representation = UnitSphericalWrap180Representation
            # Make a Spherical Wrap180 instead
            elif isinstance(self._data, SphericalRepresentation):
                self._data = SphericalWrap180Representation(
                    lat=self._data.lat, lon=self._data.lon, distance=self._data.distance)
                self.representation = SphericalWrap180Representation

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
        if isinstance(self._data, SphericalRepresentation):
            return self

        if not isinstance(self.observer, BaseCoordinateFrame):
            raise ConvertError("Cannot calculate distance to the solar disk "
                               "for observer '{}' "
                               "without `obstime` being specified.".format(self.observer))

        rep = self.represent_as(UnitSphericalWrap180Representation)
        lat, lon = rep.lat, rep.lon
        alpha = np.arccos(np.cos(lat) * np.cos(lon)).to(lat.unit)
        c = self.observer.radius**2 - self.rsun**2
        b = -2 * self.observer.radius * np.cos(alpha)
        d = ((-1*b) - np.sqrt(b**2 - 4*c)) / 2

        return self.realize_frame(SphericalWrap180Representation(lon=lon,
                                                                 lat=lat,
                                                                 distance=d))


class HelioprojectiveRadial(Helioprojective):
    """
    The Helioprojective-Radial frame is a spherical coordinate system projected
    on to the the celestial sphere.

    The center of the solar disk is defined to be at the south pole of the
    sphere and by definition has $\theta_p$ as 0 at this pole. This frame
    however, stores the declination parameter $\delta_p = $\theta_p - 90\deg$,
    so the center of the solar disk is $(0, -90)$.

    Parameters
    ----------
    dec: `~astropy.coordinates.Latitude`
        Declination Parameter.
    psi: `~astropy.coordinates.Longitude`
        Longitude coordinate.
    distance: `~astropy.units.Quantity`
        The radial distance from the observer to the coordinate point.
    """

    default_representation = SouthPoleSphericalRepresentation

    frame_specific_representation_info = {
        SouthPoleSphericalRepresentation: [RepresentationMapping(reprname='lon',
                                                                 framename='psi',
                                                                 defaultunit=u.arcsec),
                                           RepresentationMapping(reprname='lat',
                                                                 framename='el',
                                                                 defaultunit=u.arcsec),
                                           RepresentationMapping(reprname='distance',
                                                                 framename='distance',
                                                                 defaultunit=None)],

        UnitSouthPoleSphericalRepresentation: [RepresentationMapping(reprname='lon',
                                                                     framename='psi',
                                                                     defaultunit=u.arcsec),
                                               RepresentationMapping(reprname='lat',
                                                                     framename='el',
                                                                     defaultunit=u.arcsec)],

        SphericalRepresentation: [RepresentationMapping(reprname='lon',
                                                        framename='psi',
                                                        defaultunit=u.arcsec),
                                  RepresentationMapping(reprname='lat',
                                                        framename='dec',
                                                        defaultunit=u.arcsec),
                                  RepresentationMapping(reprname='distance',
                                                        framename='distance',
                                                        defaultunit=None)],

        UnitSphericalRepresentation: [RepresentationMapping(reprname='lon',
                                                            framename='psi',
                                                            defaultunit=u.arcsec),
                                      RepresentationMapping(reprname='lat',
                                                            framename='dec',
                                                            defaultunit=u.arcsec)],
    }

    obstime = TimeFrameAttributeSunPy()
    rsun = Attribute(default=RSUN_METERS.to(u.km))
    observer = ObserverCoordinateAttribute(HeliographicStonyhurst, default="earth")

    def calculate_distance(self):
        """
        This method calculates the third coordinate of the Helioprojective
        frame. It assumes that the coordinate point is on the disk of the Sun
        at the rsun radius.

        If a point in the frame is off limb then NaN will be returned.

        Returns
        -------
        new_frame : `~sunpy.coordinates.frames.HelioProjectiveRadial`
            A new frame instance with all the attributes of the original but
            now with a third coordinate.
        """
        # Skip if we already are 3D
        if isinstance(self._data, (SphericalRepresentation, SouthPoleSphericalRepresentation)):
            return self

        if not isinstance(self.observer, BaseCoordinateFrame):
            raise ConvertError("Cannot calculate distance to the solar disk "
                               "for observer '{}' "
                               "without `obstime` being specified.".format(self.observer))

        rep = self.represent_as(UnitSouthPoleSphericalRepresentation)

        distance = self.observer.radius - (self.rsun * np.cos(rep.theta))

        return self.realize_frame(SouthPoleSphericalRepresentation(phi=rep.phi,
                                                                   theta=rep.theta,
                                                                   distance=distance))
