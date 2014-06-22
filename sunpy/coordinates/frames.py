"""
SunPy's built-in coordinate frames.
Part of the proposed Coordinates API.
@author: Pritish C. (VaticanCameos)
"""

# NumPy import
import numpy as np

# Astropy imports
from astropy.extern import six
from astropy.utils.compat.odict import OrderedDict
from astropy import units as u
from astropy.coordinates.representation import (SphericalRepresentation, CylindricalRepresentation,
                                                CartesianRepresentation)
from astropy.coordinates.baseframe import BaseCoordinateFrame, frame_transform_graph
from astropy.coordinates.transformations import FunctionTransform, DynamicMatrixTransform
from astropy.coordinates import FrameAttribute

# SunPy imports
from sunpy import sun as s # For Carrington rotation number
from sunpy.sun import constants

RSUN_METERS = constants.constant('radius').si.value
DSUN_METERS = constants.constant.au.si.value

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
        the solar radius. Optional, must be keyword.
    """
    
    default_representation = SphericalRepresentation

    _frame_specific_representation_info = {
        'spherical': {'names': ('lon', 'lat', 'rad'), 'units': (u.deg, u.deg, u.km)},
        }

    rad = FrameAttribute(default=((RSUN_METERS/1000)*u.km))

def _carrington_offset():
    # This method is to return the Carrington offset.
    return s.heliographic_solar_center()[0]

class HelioGraphicCarrington(HelioGraphicStonyhurst):
    """
    A coordinate or frame in the Carrington Heliographic
    system.
    This frame differs from the Stonyhurst version in the
    definition of the longitude, which is defined using
    an offset which is a time-dependent scalar value.
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
    """
    
    default_representation = SphericalRepresentation

    _frame_specific_representation_info = {
        'spherical': {'names': ('lon', 'lat', 'rad'), 'units': (u.deg, u.deg, u.km)},
        }

    rad = FrameAttribute(default=((RSUN_METERS/1000)*u.km))

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
        X-axis coordinate, specified in kilometres. Optional, must
        be keyword.
    y: `Quantity` object.
        Y-axis coordinate, specified in kilometres. Optional, must
        be keyword.
    z: `Quantity` object. Shared by both representations.
        Z-axis coordinate, specified in kilometres. Optional, must
        be keyword.
    d: `Quantity` object.
        Represents the distance between the observer and the feature.
        Defaults to 1AU.
    D0: `Quantity` object.
        Represents the distance between the observer and the Sun center.
        Defaults to 1RSUN.
    """

    default_representation = CartesianRepresentation

    _frame_specific_representation_info = {
        'cartesian': {'names': ('x', 'y', 'z'), 'units': (u.km, u.km, u.km)},
        'cylindrical': {'names': ('rho', 'psi', 'z'), 'units': (None, u.deg, u.km)}
        }

    d = FrameAttribute(default=(1*u.au).to(u.km))
    D0 = FrameAttribute(default=((RSUN_METERS/1000)*u.km))
    
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
    zeta: Z-axis coordinate.
        Defined as D0 - d when transforming from Heliocentric
        coordinates.
    """

    default_representation = CartesianRepresentation

    _frame_specific_representation_info = {
        'cartesian': {'names': ('Tx', 'Ty', 'zeta'), 'units': (u.deg, u.deg, u.km)},
        'cylindrical': {'names': ('Trho', 'psi', 'z'), 'units': (u.deg, u.deg, u.km)}
        }

    # Note that Trho = Drho + 90, and Drho is the declination parameter.
    # According to Thompson, we use Trho internally and Drho as part of
    # the (Drho, psi) pair when defining a coordinate in this system.

# ------------------ Transformation Framework -------------------------
# This portion is reserved for the implementation of transformations
# as defined by Thompson.

@frame_transform_graph.transform(FunctionTransform, HelioGraphicStonyhurst, HelioGraphicCarrington)
def hcs_to_hcg(hcscoord, hcgframe):
    c_lon = hcscoord.spherical.lon + _carrington_offset()
    representation = SphericalRepresentation(c_lon, hcscoord.spherical.lat)
    return HelioGraphicCarrington(representation)

@frame_transform_graph.transform(FunctionTransform, HelioGraphicCarrington, HelioGraphicStonyhurst)
def hcg_to_hcs(hcgcoord, hcsframe):
    s_lon = hcgcoord.spherical.lon - _carrington_offset()
    representation = SphericalRepresentation(s_lon, hcgcoord.spherical.lat)
    return HelioGraphicStonyhurst(representation)

@frame_transform_graph.transform(FunctionTransform, HelioCentric, HelioProjective)
def helioc_to_heliop(helioccoord, heliopframe):
    # Calculate z, assuming it is on the Sun's surface.
    x = helioccoord.cartesian.x.value * 1000
    y = helioccoord.cartesian.y.value * 1000
    z = helioccoord.cartesian.z.value * 1000
    zeta = DSUN_METERS - z

    distance = np.sqrt(x ** 2 + y ** 2 + zeta ** 2)
    hpcx = np.rad2deg(np.arctan2(x, zeta))
    hpcy = np.rad2deg(np.arcsin(y / distance))

    representation = CartesianRepresentation(hpcx, hpcy, zeta)
    return HelioProjective(representation)
    
@frame_transform_graph.transform(FunctionTransform, HelioProjective, HelioCentric)
def heliop_to_helioc(heliopcoord, heliocframe):
    x = heliopcoord.cartesian.x.value
    y = heliopcoord.cartesian.y.value
    c = np.array([np.deg2rad(1), np.deg2rad(1)])

    cosx = np.cos(x * c[0])
    sinx = np.sin(x * c[0])
    cosy = np.cos(y * c[1])
    siny = np.sin(y * c[1])
    
    q = DSUN_METERS * cosy * cosx
    distance = q ** 2 - DSUN_METERS ** 2 + RSUN_METERS ** 2
    distance = q - np.sqrt(distance)

    rx = distance * cosy * sinx
    ry = distance * siny
    rz = DSUN_METERS - distance * cosy * cosx

    representation = CartesianRepresentation(rx, ry, rz)
    return HelioCentric(representation)

@frame_transform_graph.transform(FunctionTransform, HelioCentric, HelioGraphicStonyhurst)
def hcc_to_hgs(helioccoord, heliogframe):
    x = helioccoord.cartesian.x.value * 1000
    y = helioccoord.cartesian.y.value * 1000
    z = helioccoord.cartesian.z.value * 1000
    
    l0_deg = _carrington_offset()
    b0_deg = s.heliographic_solar_center()[1]

    cosb = np.cos(np.deg2rad(b0_deg))
    sinb = np.sin(np.deg2rad(b0_deg))

    hecr = np.sqrt(x**2 + y**2 + z**2)
    hgln = np.arctan2(x, z * cosb - y * sinb) + np.deg2rad(l0_deg)
    hglt = np.arcsin((y * cosb + z * sinb) / hecr)

    representation = SphericalRepresentation(np.rad2deg(hgln),
                                             np.rad2deg(hglt),
                                             hecr)
    return HelioGraphicStonyhurst(representation)

@frame_transform_graph.transfor(FunctionTransform, HelioGraphicStonyhurst, HelioCentric)
def hgs_to_hcc(heliogcoord, heliopframe):
    hglon = heliogcoord.spherical.lon.value
    hglat = heliogcoord.spherical.lat.value
    r = heliogcoord.spherical.distance.value * 1000

    l0_deg = _carrington_offset()
    b0_deg = s.heliographic_solar_center()[1]

    lon = np.deg2rad(hglon)
    lat = np.deg2rad(hglat)

    cosb = np.cos(np.deg2rad(b0_deg))
    sinb = np.sin(np.deg2rad(b0_deg))

    lon = lon - np.deg2rad(l0_deg)

    cosx = np.cos(lon)
    sinx = np.sin(lon)
    cosy = np.cos(lat)
    siny = np.sin(lat)

    x = r * cosy * sinx
    y = r * (siny * cosb - cosy * cosx * sinb)
    zz = r * (siny * sinb + cosy * cosx * cosb)

    representation = CartesianRepresentation(x, y, zz)
    return HelioCentric(representation)

