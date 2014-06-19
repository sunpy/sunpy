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

# SunPy imports
from sunpy import sun as s # For Carrington rotation number
from sunpy.sun import constants

RADIUS = constants.constant('radius')

class HelioGraphicStonyhurst(BaseCoordinateFrame):
    """
    A coordinate or frame in the Stonyhurst Heliographic
    system.
    This system is known to remain fixed with respect to
    the center of the Earth, and its quantities, the
    latitude and longitude, are specified in degrees.

    Parameters
    ----------
    representation: `BaseRepresentation` or None
        A representation object or None to have no data.
    lon: `Angle` object.
        The longitude for this object (``lat`` must also be given and ``representation``
        must be None).
    lat: `Angle` object.
        The latitude for this object (``lon`` must also be given and ``representation``
        must be None).
    rad: `astropy.units.Quantity` object.
        This quantity holds the radial distance. If not specified, it is, by default,
        the solar radius.
    """
    default_representation = SphericalRepresentation

    frame_attr_names = {}

    _frame_specific_representation_info = {
        'spherical': {'names': ('lon', 'lat', 'rad'), 'units': (u.deg, u.deg, u.km)},
        }

    def __init__(self, *args, **kwargs):
        if kwargs is None:
            # If only lon and lat are specified.
            if len(args) < 3:
                args.append(RADIUS*u.km)
        elif args is None:
            if 'rad' not in kwargs:
                kwargs['rad'] = RADIUS*u.km
        else:
            # Any other cases?
            pass
        super(HelioGraphicStonyhurst, self).__init__(args, kwargs)

def _carrington_offset():
    # This method is to return the Carrington offset.
    pass

class HelioGraphicCarrington(HelioGraphicStonyhurst):
    """
    A coordinate or frame in the Carrington Heliographic
    system.
    This frame differs from the Stonyhurst version in the
    definition of the longitude, which is defined using
    an offset which is a time-dependent scalar value.
    representation: `BaseRepresentation` or None.
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
        the solar radius.
    """
    default_representation = SphericalRepresentation

    frame_attr_names = {}

    _frame_specific_representation_info = {
        'spherical': {'names': ('lon', 'lat', 'rad'), 'units': (u.deg, u.deg, u.km)},
        }

    def __init__(self, *args, **kwargs):

        # TODO: Understand how Carrington offset works.
        super(HelioGraphicCarrington, self).__init__(args, kwargs)

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
    representation: `BaseRepresentation` or None.
        A representation object. If specified, other parameters must
        be in keyword form.
    x: `Quantity` object.
        X-axis coordinate, specified in kilometres.
    y: `Quantity` object.
        Y-axis coordinate, specified in kilometres.
    z: `Quantity` object. Shared by both representations.
        Z-axis coordinate, specified in kilometres.
    """

    default_representation = CartesianRepresentation

    frame_attr_names = {}

    _frame_specific_representation_info = {
        'cartesian': {'names': ('x', 'y', 'z'), 'units': (u.km, u.km, u.km)},
        'cylindrical': {'names': ('rho', 'psi', 'z'), 'units': (None, u.deg, u.km)}
        }
    
