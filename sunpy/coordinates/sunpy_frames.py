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

RSUN = constants.constant('radius')

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

    rad = FrameAttribute(default=((RSUN.value/1000)*u.km))

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

    rad = FrameAttribute(default=((RSUN.value/1000)*u.km))

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
    D0 = FrameAttribute(default=((RSUN.value/1000)*u.km))
    
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
        'cylindrical': {'names': ('Trho', 'psi', 'z'), 'units': (u.deg, u.deg, None)}
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
    mult_factor = 180/(np.pi * helioccoord.d.value)
    Tx = mult_factor * helioccoord.cartesian.x.value * u.deg
    Ty = mult_factor * helioccoord.cartesian.y.value * u.deg
    zeta = helioccoord.D0 - helioccoord.d
    representation = CartesianRepresentation(Tx, Ty, zeta)
    return HelioProjective(representation)

@frame_transform_graph.transform(FunctionTransform, HelioProjective, HelioCentric)
def heliop_to_helioc(heliopcoord, heliocframe):
    # mult_factor = (np.pi * heliocframe.d.value)/180
    # x = mult_factor * heliopcoord.cartesian.x.value * u.km
    # y = mult_factor * heliopcoord.cartesian.y.value * u.km
    pass
    
