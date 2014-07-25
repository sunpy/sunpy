"""
SunPy's built-in coordinate frames.
Part of the proposed Coordinates API.
@author: Pritish C. (VaticanCameos)
"""

# NumPy import
import numpy as np

# Astropy imports
from astropy import units as u
from astropy.coordinates.representation import (BaseRepresentation,
                                                CartesianRepresentation)
from astropy.coordinates.baseframe import (BaseCoordinateFrame, frame_transform_graph,
                                           RepresentationMapping)
from astropy.coordinates.transformations import FunctionTransform
from astropy.coordinates import FrameAttribute

# SunPy imports
from sunpy import sun as s # For Carrington rotation number
from representation import SphericalWrap180Representation

from datetime import datetime

from frameattributes import TimeFrameAttributeSunPy

RSUN_METERS = s.constants.constant('radius').si
DSUN_METERS = s.constants.constant('mean distance').si

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
    hlon: `Angle` object.
        The longitude for this object (``lat`` must also be given and ``representation``
        must be None).
    hlat: `Angle` object.
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
    hlon=1.0 deg, hlat=1.0 deg, rad=2.0 km>
    >>> sc.frame
    <HelioGraphicStonyhurst Coordinate: dateobs=2010-01-01 00:00:45, 
    hlon=1.0 deg, hlat=1.0 deg, rad=2.0 km>
    >>> sc = SkyCoord(HelioGraphicStonyhurst(-10*u.deg, 2*u.deg))
    >>> sc
    <SkyCoord (HelioGraphicStonyhurst): dateobs=None, hlon=-10.0 deg, 
    hlat=2.0 deg, rad=695508.0 km>
    """

    default_representation = SphericalWrap180Representation

    _frame_specific_representation_info = {
        'sphericalwrap180': [RepresentationMapping('lon', 'hlon', 'recommended'),
                             RepresentationMapping('lat', 'hlat', 'recommended'),
                             RepresentationMapping('distance', 'rad', 'recommended')],
        }

    dateobs = TimeFrameAttributeSunPy()

    def __init__(self, *args, **kwargs):
        if args or kwargs: # Non-empty frame use case.
            if args and kwargs: # Mixed use case.
                if args and not isinstance(args[0], BaseRepresentation):
                    if len(args) > 0 and len(args) <= 2 and 'rad' not in kwargs:
                    # If one of hlon/hlat are in args
                        kwargs['rad'] = kwargs.get('rad', RSUN_METERS.to(u.km))
            elif not args: # kwargs-only use case
                if 'representation' not in kwargs:
                    #if 'rad' not in kwargs: # This default is required by definition.
                    if 'hlon' in kwargs and 'hlat' in kwargs:
                        kwargs['rad'] = kwargs.get('rad', RSUN_METERS.to(u.km))
            elif not kwargs: # args-only use case.
                if len(args) == 2:
                    args = list(args)
                    args.append(RSUN_METERS.to(u.km))
                    args = tuple(args)

        super(HelioGraphicStonyhurst, self).__init__(*args, **kwargs)

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
    hlon: `Angle` object.
        The longitude for this object (``lat`` must also be given and ``representation``
        must be None).
    hlat: `Angle` object.
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
    hlon=1.0 deg, hlat=2.0 deg, rad=3.0 km>
    >>> sc = SkyCoord([1,2,3]*u.deg, [4,5,6]*u.deg, [5,6,7]*u.km, 
    dateobs="2010/01/01T00:00:45", frame="heliographiccarrington")
    >>> sc
    <SkyCoord (HelioGraphicCarrington): dateobs=2010-01-01 00:00:45, 
    (hlon, hlat, rad) in (deg, deg, km)
        [(1.0, 4.0, 5.0), (2.0, 5.0, 6.0), (3.0, 6.0, 7.0)]>
    """

    default_representation = SphericalWrap180Representation

    _frame_specific_representation_info = {
        'sphericalwrap180': [RepresentationMapping('lon', 'hlon', 'recommended'),
                             RepresentationMapping('lat', 'hlat', 'recommended'),
                             RepresentationMapping('distance', 'rad', 'recommended')]
        }

    #rad = FrameAttribute(default=((RSUN_METERS/1000)*u.km))
    dateobs = TimeFrameAttributeSunPy()

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
        X-axis coordinate, specified in degrees.
    Ty: `Angle` object.
        Y-axis coordinate, specified in degrees.
    distance: Z-axis coordinate.
        Represents the radial distance between the solar center
        and the observer.
        Defaults to 1AU.
    zeta: `Quantity` object.
        Represents the distance between observer and feature/point.
        Defaults to 0.
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
    >>> hp = HelioProjective(0*u.deg, 0*u.deg, zeta=1*u.km, 
    dateobs="2010/01/01T00:00:00")
    >>> hp
    <HelioProjective Coordinate: dateobs=2010-01-01 00:00:00, D0=149597870.7 km
    , Tx=0.0 arcsec, Ty=0.0 arcsec, distance=149597869.7 km>
    >>> sc = SkyCoord(hp)
    >>> sc
    <SkyCoord (HelioProjective): dateobs=2010-01-01 00:00:00, D0=149597870.7 km
    , Tx=0.0 arcsec, Ty=0.0 arcsec, distance=149597869.7 km>  
    """

    default_representation = SphericalWrap180Representation

    _frame_specific_representation_info = {
        'sphericalwrap180': [RepresentationMapping('lon', 'Tx', u.arcsec),
                             RepresentationMapping('lat', 'Ty', u.arcsec),
                             RepresentationMapping('distance', 'distance', u.km)],
        'cylindrical': [RepresentationMapping('rho', 'Trho', u.arcsec),
                        RepresentationMapping('phi', 'psi', u.arcsec),
                        RepresentationMapping('distance', 'distance', u.km)]}

    D0 = FrameAttribute(default=(1*u.au).to(u.km))
    dateobs = TimeFrameAttributeSunPy()

    @property
    def zeta(self):
        """zeta is defined as a property."""
        return self.D0 - self.distance

    def __init__(self, *args, **kwargs):
        """
        This is the custom constructor method for HelioProjective frames.
        It is required as we wish to default 'distance' to D0 - zeta when it
        itself is not present and there are no supporting arguments.
        'zeta' is a supporting argument that must be specified as a kwarg.
        If 'zeta' is present, 'distance' can be calculated as given.
        Both 'zeta' and 'distance' cannot be present at the same time.
        """
        if args or (kwargs and len(kwargs) != 1):
            # Non-empty frame use case.
            if args and kwargs:
                # If we have both args and kwargs.
                if isinstance(args[0], BaseRepresentation):
                    # The case when first arg is a representation.
                    if 'zeta' in kwargs:
                        # zeta cannot be provided as SphericalRep takes three arguments.
                        raise TypeError("zeta cannot be specified with a representation "
                                        "for the {0} frame.".format(self.__class__))
                elif len(args) < 3:
                    # If we have either args(Tx) and rest kwargs, or args(Tx, Ty) and rest kwargs.
                    if 'distance' not in kwargs and 'zeta' in kwargs:
                        kwargs['distance'] = kwargs.get('D0', self.D0) - kwargs['zeta']
                        kwargs.pop('zeta')
                    elif 'distance' not in kwargs and 'zeta' not in kwargs:
                        kwargs['distance'] = self.get_distance_hpc(*args, **kwargs)
                    elif 'distance' in kwargs and 'zeta' in kwargs:
                        raise TypeError("zeta and distance cannot both be "
                                        "specified in the {0} frame.".format(self.__class__.name))
                elif len(args) == 3:
                    # If we have args(Tx, Ty, distance).
                    if 'zeta' in kwargs:
                        raise TypeError("zeta and distance cannot both "
                                        "be specified here for the {0} frame.".format(self.__class__.name))
            elif not kwargs:
                # The case when kwargs are not present.
                if len(args) == 2:
                    # args(Tx, Ty) provided.
                    args = list(args)
                    args.append(self.get_distance_hpc(*tuple(args)))
                    args = tuple(args)
            elif not args:
                # The case when args are not present.
                if 'distance' not in kwargs and 'zeta' in kwargs:
                    kwargs['distance'] = kwargs.get('D0', self.D0) - kwargs['zeta']
                    kwargs.pop('zeta')
                elif 'distance' not in kwargs and 'zeta' not in kwargs:
                    if 'Tx' in kwargs and 'Ty' in kwargs or 'lon' in kwargs and 'lat' in kwargs:
                        # This if clause was added to deal with a frame
                        # which does not have Tx, Ty, distance but may
                        # have other kwargs (FrameAttributes).
                        kwargs['distance'] = self.get_distance_hpc(**kwargs)
                elif 'distance' in kwargs and 'zeta' in kwargs:
                    raise TypeError("zeta and distance cannot both be "
                                    "specified here for the {0} frame.".format(self.__class__))
        # Finally, make a call to the super constructor.
        super(HelioProjective, self).__init__(*args, **kwargs)
           
    # Note that Trho = Drho + 90, and Drho is the declination parameter.
    # According to Thompson, we use Trho internally and Drho as part of
    # the (Drho, psi) pair when defining a coordinate in this system.
           
    def get_distance_hpc(self, *args, **kwargs):
        """
        This method calculates the 'distance' parameter if it is not specifed.
        It takes Tx, Ty and possibly D0 as its input.
        It returns distance in kilometers.
        """
        c, b = None, None
        list_params = self._get_input_params(*args, **kwargs)        
        
        alpha = np.arccos(np.cos(list_params[0]) * np.cos(list_params[1]))\
                .to(list_params[0].unit)
        if len(list_params) == 3:
            c = (list_params[2].to(u.m))**2 - RSUN_METERS**2
            b = -2 * list_params[2].to(u.m) * np.cos(alpha)
        else:
            c = (self.D0.to(u.m))**2 - RSUN_METERS**2
            b = -2 * self.D0.to(u.m) * np.cos(alpha)
        d = ((-1*b) - np.sqrt(b**2 - 4*c)) / 2
        
        return d.to(u.km)

    def _get_input_params(self, *args, **kwargs):
        """
        This method prunes the args and kwargs to find Tx, Ty and D0.
        It is used by get_distance_hpc(). It returns a list of these values
        in the order [Tx, Ty, D0] or [Tx, Ty].
        If D0 is not specified, get_distance_hpc() uses the default D0 instead.
        """
        list_params = []
        if len(args) == 2:
            list_params = list(args)
        elif len(args) == 1:
            if 'Tx' in kwargs:
                list_params.append(kwargs['Tx'])
                list_params.append(args[0])
            elif 'Ty' in kwargs:
                list_params = list(args)
                list_params.append(kwargs['Ty'])
        elif not args:
            list_params.append(kwargs['Tx'])
            list_params.append(kwargs['Ty'])
        if 'D0' in kwargs:
            list_params.append(kwargs['D0'])
        return list_params
        
def _carrington_offset(dateobs):
    if dateobs is None:
        raise ValueError("To perform this transformation the coordinate Frame needs a dateobs Attribute")
    # This method is to return the Carrington offset.
    return s.heliographic_solar_center(dateobs)[0]


# ------------------ Transformation Framework -------------------------
# This portion is reserved for the implementation of transformations
# as defined by Thompson.

@frame_transform_graph.transform(FunctionTransform, HelioGraphicStonyhurst, HelioGraphicCarrington)
def hgs_to_hgc(hgscoord, hgcframe):
    c_lon = hgscoord.spherical.lon + _carrington_offset(hgscoord.dateobs).to(u.deg)
    representation = SphericalWrap180Representation(c_lon, hgscoord.hlat, hgscoord.rad)
    return hgcframe.realize_frame(representation)

@frame_transform_graph.transform(FunctionTransform, HelioGraphicCarrington, HelioGraphicStonyhurst)
def hgc_to_hgs(hgccoord, hgsframe):
    s_lon = hgccoord.spherical.lon - _carrington_offset(hgccoord.dateobs).to(u.deg)
    representation = SphericalWrap180Representation(s_lon, hgccoord.hlat, hgccoord.rad)
    return hgsframe.realize_frame(representation)

@frame_transform_graph.transform(FunctionTransform, HelioCentric, HelioProjective)
def helioc_to_heliop(helioccoord, heliopframe):
    x = helioccoord.x.to(u.m)
    y = helioccoord.y.to(u.m)
    z = helioccoord.z.to(u.m)

    # d is calculated as the distance between the points
    # (x,y,z) and (0,0,D0).
    distance = np.sqrt(x**2 + y**2 + (helioccoord.D0.to(u.m) - z)**2)

    hpcx = np.rad2deg(np.arctan2(x, helioccoord.D0 - z))
    hpcy = np.rad2deg(np.arcsin(y / distance))

    representation = SphericalWrap180Representation(hpcx, hpcy, distance.to(u.km))
    return heliopframe.realize_frame(representation)

@frame_transform_graph.transform(FunctionTransform, HelioProjective, HelioCentric)
def heliop_to_helioc(heliopcoord, heliocframe):
    x = np.deg2rad(heliopcoord.Tx)
    y = np.deg2rad(heliopcoord.Ty)

    cosx = np.cos(x)
    sinx = np.sin(x)
    cosy = np.cos(y)
    siny = np.sin(y)

    rx = (heliopcoord.distance.to(u.m)) * cosy * sinx
    ry = (heliopcoord.distance.to(u.m)) * siny
    rz = (heliopcoord.D0.to(u.m)) - (heliopcoord.distance.to(u.m)) * cosy * cosx

    representation = CartesianRepresentation(rx.to(u.km), ry.to(u.km), rz.to(u.km))
    return heliocframe.realize_frame(representation)

@frame_transform_graph.transform(FunctionTransform, HelioCentric, HelioGraphicStonyhurst)
def hcc_to_hgs(helioccoord, heliogframe):
    x = helioccoord.x.to(u.m)
    y = helioccoord.y.to(u.m)
    z = helioccoord.z.to(u.m)

    l0b0_pair = s.heliographic_solar_center()

    l0_rad = l0b0_pair[0]
    b0_deg = l0b0_pair[1]

    cosb = np.cos(np.deg2rad(b0_deg))
    sinb = np.sin(np.deg2rad(b0_deg))

    hecr = np.sqrt(x**2 + y**2 + z**2)
    hgln = np.arctan2(x, z * cosb - y * sinb) + l0_rad
    hglt = np.arcsin((y * cosb + z * sinb) / hecr)

    representation = SphericalWrap180Representation(np.rad2deg(hgln),
                                             np.rad2deg(hglt),
                                             hecr.to(u.km))
    return heliogframe.realize_frame(representation)

@frame_transform_graph.transform(FunctionTransform, HelioGraphicStonyhurst, HelioCentric)
def hgs_to_hcc(heliogcoord, heliocframe):
    hglon = heliogcoord.hlon
    hglat = heliogcoord.hlat
    r = heliogcoord.rad.to(u.m)
    
    l0b0_pair = s.heliographic_solar_center()

    l0_rad = l0b0_pair[0]
    b0_deg = l0b0_pair[1]

    lon = np.deg2rad(hglon)
    lat = np.deg2rad(hglat)

    cosb = np.cos(b0_deg.to(u.rad))
    sinb = np.sin(b0_deg.to(u.rad))

    lon = lon - l0_rad

    cosx = np.cos(lon)
    sinx = np.sin(lon)
    cosy = np.cos(lat)
    siny = np.sin(lat)

    x = r * cosy * sinx
    y = r * (siny * cosb - cosy * cosx * sinb)
    zz = r * (siny * sinb + cosy * cosx * cosb)

    representation = CartesianRepresentation(x.to(u.km), y.to(u.km), zz.to(u.km))
    return heliocframe.realize_frame(representation)
