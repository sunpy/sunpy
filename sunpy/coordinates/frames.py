"""
SunPy's built-in coordinate frames.
Part of the proposed Coordinates API.
@author: Pritish C. (VaticanCameos)
"""

# NumPy import
import numpy as np

# Astropy imports
from astropy import units as u
from astropy.coordinates.representation import (CartesianRepresentation,
                                                UnitSphericalRepresentation,
                                                SphericalRepresentation)
from astropy.coordinates.baseframe import (BaseCoordinateFrame,
                                           RepresentationMapping)
from astropy.coordinates import FrameAttribute

# SunPy imports
from sunpy import sun # For Carrington rotation number
from representation import (SphericalWrap180Representation,
                            UnitSphericalWrap180Representation)

from .frameattributes import TimeFrameAttributeSunPy

RSUN_METERS = sun.constants.constant('radius').si.to(u.m)
DSUN_METERS = sun.constants.constant('mean distance').si.to(u.m)

__all__ = ['HelioGraphicStonyhurst', 'HelioGraphicCarrington',
           'HelioCentric', 'HelioProjective']

class HelioGraphicStonyhurst(BaseCoordinateFrame):
    """
    A coordinate or frame in the Stonyhurst Heliographic
    system.
    This system is known to remain fixed with respect to
    the center of the Earth, and its quantities, the
    latitude and longitude, are specified in degrees.

    Parameters
    ----------
    representation: `~astropy.coordinates.BaseRepresentation` or None
        A representation object or None to have no data.
    lon: `Angle` object.
        The longitude for this object (``lat`` must also be given and ``representation``
        must be None).
    lat: `Angle` object.
        The latitude for this object (``lon`` must also be given and ``representation``
        must be None).
    rad: `astropy.units.Quantity` object.
        This quantity holds the radial distance. If not specified, it is, by default,
        the solar radius. Optional, must be keyword

    Examples
    --------
    >>> sc = SkyCoord(1*u.deg, 1*u.deg, 2*u.km, frame="heliographicstonyhurst",
    dateobs="2010/01/01T00:00:45")
    >>> sc
    <SkyCoord (HelioGraphicStonyhurst): dateobs=2010-01-01 00:00:45,
    lon=1.0 deg, lat=1.0 deg, rad=2.0 km>
    >>> sc.frame
    <HelioGraphicStonyhurst Coordinate: dateobs=2010-01-01 00:00:45,
    lon=1.0 deg, lat=1.0 deg, rad=2.0 km>
    >>> sc = SkyCoord(HelioGraphicStonyhurst(-10*u.deg, 2*u.deg))
    >>> sc
    <SkyCoord (HelioGraphicStonyhurst): dateobs=None, lon=-10.0 deg,
    lat=2.0 deg, rad=695508.0 km>

    Notes
    -----
    This frame will always be converted a 3D frame where the radius defaults to
    RSun.
    """

    default_representation = SphericalWrap180Representation

    _frame_specific_representation_info = {
        'spherical': [RepresentationMapping('lon', 'lon', 'recommended'),
                      RepresentationMapping('lat', 'lat', 'recommended'),
                      RepresentationMapping('distance', 'rad', 'recommended')],

        'sphericalwrap180': [RepresentationMapping('lon', 'lon', 'recommended'),
                             RepresentationMapping('lat', 'lat', 'recommended'),
                             RepresentationMapping('distance', 'rad', 'recommended')]
        }

    dateobs = TimeFrameAttributeSunPy()
    RSun = FrameAttribute(default=RSUN_METERS.to(u.km))

    def __init__(self, *args, **kwargs):
        _rep_kwarg = kwargs.get('representation', None)

        super(HelioGraphicStonyhurst, self).__init__(*args, **kwargs)

        # If representation was explicitly passed, do not change the rep.
        if not _rep_kwarg:
            # The base __init__ will make this a UnitSphericalRepresentation
            # This makes it Wrap180 instead
            if isinstance(self._data, UnitSphericalRepresentation):
                self._data = SphericalWrap180Representation(lat=self._data.lat,
                                                            lon=self._data.lon,
                                                            distance=self.RSun)
                self.representation = SphericalWrap180Representation

            # Make a Spherical Wrap180 instead
            elif isinstance(self._data, SphericalRepresentation):
                self._data = SphericalWrap180Representation(lat=self._data.lat,
                                                            lon=self._data.lon,
                                                            distance=self._data.distance)
                self.representation = SphericalWrap180Representation




class HelioGraphicCarrington(HelioGraphicStonyhurst):
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
        The longitude for this object (``lat`` must also be given and ``representation``
        must be None).
    lat: `Angle` object.
        The latitude for this object (``lon`` must also be given and ``representation``
        must be None).
    rad: `astropy.units.Quantity` object, optional, must be keyword.
        This quantity holds the radial distance. If not specified, it is, by default,
        the solar radius. Optional, must be keyword.

    Examples
    --------
    >>> sc = SkyCoord(1*u.deg, 2*u.deg, 3*u.km, frame="heliographiccarrington",
    dateobs="2010/01/01T00:00:30")
    >>> sc
    <SkyCoord (HelioGraphicCarrington): dateobs=2010-01-01 00:00:30,
    lon=1.0 deg, lat=2.0 deg, rad=3.0 km>
    >>> sc = SkyCoord([1,2,3]*u.deg, [4,5,6]*u.deg, [5,6,7]*u.km,
    dateobs="2010/01/01T00:00:45", frame="heliographiccarrington")
    >>> sc
    <SkyCoord (HelioGraphicCarrington): dateobs=2010-01-01 00:00:45,
    (lon, lat, rad) in (deg, deg, km)
        [(1.0, 4.0, 5.0), (2.0, 5.0, 6.0), (3.0, 6.0, 7.0)]>
    """

    default_representation = SphericalWrap180Representation

    _frame_specific_representation_info = {
        'spherical': [RepresentationMapping('lon', 'lon', 'recommended'),
                      RepresentationMapping('lat', 'lat', 'recommended'),
                      RepresentationMapping('distance', 'rad', 'recommended')],

        'sphericalwrap180': [RepresentationMapping('lon', 'lon', 'recommended'),
                             RepresentationMapping('lat', 'lat', 'recommended'),
                             RepresentationMapping('distance', 'rad', 'recommended')]
        }

    dateobs = TimeFrameAttributeSunPy()
    RSun = FrameAttribute(default=RSUN_METERS.to(u.km))

class HelioCentric(BaseCoordinateFrame):
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
    D0: `Quantity` object.
        Represents the distance between the observer and the Sun center.
        Defaults to 1AU.

    Examples
    --------
    >>> sc = SkyCoord(CartesianRepresentation(10*u.km, 1*u.km, 2*u.km),
    dateobs="2011/01/05T00:00:50", frame="heliocentric")
    >>> sc
    <SkyCoord (HelioCentric): dateobs=2011-01-05 00:00:50, D0=149597870.7 km,
    x=10.0 km, y=1.0 km, z=2.0 km>
    >>> sc = SkyCoord([1,2]*u.km, [3,4]*u.m, [5,6]*u.cm, frame="heliocentric",
    dateobs="2011/01/01T00:00:54")
    >>> sc
    <SkyCoord (HelioCentric): dateobs=2011-01-01 00:00:54, D0=149597870.7 km,
    (x, y, z) in (km, m, cm)
        [(1.0, 3.0, 5.0), (2.0, 4.0, 6.0)]>
    """

    default_representation = CartesianRepresentation

    _frame_specific_representation_info = {
        'cylindrical': [RepresentationMapping('phi', 'psi', u.deg)]}

   # d = FrameAttribute(default=(1*u.au).to(u.km))
    D0 = FrameAttribute(default=(1*u.au).to(u.km))
    dateobs = TimeFrameAttributeSunPy()
    L0 = FrameAttribute(default=0*u.deg)
    B0 = FrameAttribute(default=0*u.deg)

class HelioProjective(BaseCoordinateFrame):
    """
    A coordinate or frame in the Helioprojective
    system.
    This is the projected equivalent of the Heliocentric
    coordinate system. As such, the Cartesian representation
    has degrees for each of the units, and the cylindrical
    representation has the rho parameter replaced by Trho,
    or theta_rho.

    Parameters
    ----------
    representation: `~astropy.coordinates.BaseRepresentation` or None.
        A representation object. If specified, other parameters must
        be in keyword form.
    Tx: `Angle` object.
        X-axis coordinate.
    Ty: `Angle` object.
        Y-axis coordinate.
    distance: Z-axis coordinate.
        The radial distance from the observer to the coordinate point.
    D0: `Quantity` object.
        Represents the distance between observer and solar center.
        Defaults to 1AU.

    Examples
    --------
    >>> sc = SkyCoord(0*u.deg, 0*u.deg, 5*u.km, dateobs="2010/01/01T00:00:00",
    frame="helioprojective")
    >>> sc
    <SkyCoord (HelioProjective): dateobs=2010-01-01 00:00:00, D0=149597870.7 km
    , Tx=0.0 arcsec, Ty=0.0 arcsec, distance=5.0 km>
    >>> sc = SkyCoord(0*u.deg, 0*u.deg, dateobs="2010/01/01T00:00:00",
    frame="helioprojective")
    >>> sc
    <SkyCoord (HelioProjective): dateobs=2010-01-01 00:00:00, D0=149597870.7 km
    , Tx=0.0 arcsec, Ty=0.0 arcsec, distance=149597870.7 km>
    """

    default_representation = SphericalWrap180Representation

    _frame_specific_representation_info = {
        'spherical': [RepresentationMapping('lon', 'Tx', u.arcsec),
                      RepresentationMapping('lat', 'Ty', u.arcsec),
                      RepresentationMapping('distance', 'distance', u.km)],

        'sphericalwrap180': [RepresentationMapping('lon', 'Tx', u.arcsec),
                             RepresentationMapping('lat', 'Ty', u.arcsec),
                             RepresentationMapping('distance', 'distance', u.km)],

        'unitspherical': [RepresentationMapping('lon', 'Tx', u.arcsec),
                          RepresentationMapping('lat', 'Ty', u.arcsec)],

        'unitsphericalwrap180': [RepresentationMapping('lon', 'Tx', u.arcsec),
                                 RepresentationMapping('lat', 'Ty', u.arcsec)],

        'cylindrical': [RepresentationMapping('rho', 'rho', u.km),
                        RepresentationMapping('phi', 'psi', u.arcsec),
                        RepresentationMapping('z', 'distance', u.km)]}

    D0 = FrameAttribute(default=(1*u.au).to(u.km))
    dateobs = TimeFrameAttributeSunPy()
    L0 = FrameAttribute(default=0*u.deg)
    B0 = FrameAttribute(default=0*u.deg)
    RSun = FrameAttribute(default=RSUN_METERS.to(u.km))

    def __init__(self, *args, **kwargs):
        _rep_kwarg = kwargs.get('representation', None)

        BaseCoordinateFrame.__init__(self, *args, **kwargs)

        # If representation was explicitly passed, do not change the rep.
        if not _rep_kwarg:
            # The base __init__ will make this a UnitSphericalRepresentation
            # This makes it Wrap180 instead
            if isinstance(self._data, UnitSphericalRepresentation):
                self._data = UnitSphericalWrap180Representation(lat=self._data.lat,
                                                                lon=self._data.lon)
                self.representation = UnitSphericalWrap180Representation
            # Make a Spherical Wrap180 instead
            elif isinstance(self._data, SphericalRepresentation):
                self._data = SphericalWrap180Representation(lat=self._data.lat,
                                                            lon=self._data.lon,
                                                            distance=self._data.distance)
                self.representation = SphericalWrap180Representation

    # Note that Trho = Drho + 90, and Drho is the declination parameter.
    # According to Thompson, we use Trho internally and Drho as part of
    # the (Drho, psi) pair when defining a coordinate in this system.
    def calculate_distance(self):
        """
        This method calculates the third coordnate of the Helioprojective frame.
        It assumes that the coordinate point is on the disk of the Sun at the
        RSun radius.

        If a point in the frame is off limb then NaN will be returned.

        Returns
        -------
        new_frame : `~sunpy.coordinates.frames.HelioProjective`
            A new frame instance with all the attributes of the original but now
            with a third coordinate.
        """
        # Skip if we already are 3D
        if isinstance(self._data, SphericalRepresentation):
            return self

        rep = self.represent_as(UnitSphericalWrap180Representation)
        lat, lon = rep.lat, rep.lon
        alpha = np.arccos(np.cos(lat) * np.cos(lon)).to(lat.unit)
        c = self.D0**2 - self.RSun**2
        b = -2 * self.D0.to(u.m) * np.cos(alpha)
        d = ((-1*b) - np.sqrt(b**2 - 4*c)) / 2
        return self.realize_frame(SphericalWrap180Representation(lon=lon,
                                                                 lat=lat,
                                                                 distance=d))
