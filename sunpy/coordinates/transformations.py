# -*- coding: utf-8 -*-
"""
Coordinate Transformation Functions

This module contains the functions for converting one
`sunpy.coordinates.frames` object to another.

.. warning::

  The functions in this submodule should never be called directly, transforming
  between coordinate frames should be done using the ``.transform_to`` methods
  on `~astropy.coordinates.BaseCoordinateFrame` or
  `~astropy.coordinates.SkyCoord` instances.

"""
from __future__ import absolute_import, division

import numpy as np

from astropy import units as u
from astropy.coordinates.representation import (CartesianRepresentation,
                                                UnitSphericalRepresentation,
                                                SphericalRepresentation)
from astropy.coordinates.baseframe import frame_transform_graph
from astropy.coordinates.builtin_frames import _make_transform_graph_docs
from astropy.coordinates.transformations import FunctionTransform, DynamicMatrixTransform
from astropy.coordinates.matrix_utilities import rotation_matrix, matrix_product, matrix_transpose
from astropy.coordinates import HCRS, get_body_barycentric, BaseCoordinateFrame, ConvertError

from sunpy import sun

from .representation import (SphericalWrap180Representation,
                             UnitSphericalWrap180Representation)
from .frames import (HeliographicStonyhurst, HeliographicCarrington,
                     Heliocentric, Helioprojective)

__all__ = ['hgs_to_hgc', 'hgc_to_hgs', 'hcc_to_hpc',
           'hpc_to_hcc', 'hcc_to_hgs', 'hgs_to_hcc',
           'hpc_to_hpc',
           'hcrs_to_hgs', 'hgs_to_hcrs']


def _carrington_offset(obstime):
    """
    Calculate the HG Longitude offest based on a time
    """
    if obstime is None:
        raise ValueError("To perform this transformation the coordinate"
                         " Frame needs a obstime Attribute")

    # Import here to avoid a circular import
    from .ephemeris import get_sun_L0
    return get_sun_L0(obstime)

# =============================================================================
# ------------------------- Transformation Framework --------------------------
# =============================================================================


@frame_transform_graph.transform(FunctionTransform, HeliographicStonyhurst,
                                 HeliographicCarrington)
def hgs_to_hgc(hgscoord, hgcframe):
    """
    Transform from Heliographic Stonyhurst to Heliograpic Carrington.
    """
    c_lon = hgscoord.spherical.lon + _carrington_offset(hgscoord.obstime).to(
        u.deg)
    representation = SphericalRepresentation(c_lon, hgscoord.lat,
                                                    hgscoord.radius)
    hgcframe = hgcframe.__class__(obstime=hgscoord.obstime)
    return hgcframe.realize_frame(representation)


@frame_transform_graph.transform(FunctionTransform, HeliographicCarrington,
                                 HeliographicStonyhurst)
def hgc_to_hgs(hgccoord, hgsframe):
    """
    Convert from Heliograpic Carrington to Heliographic Stonyhurst.
    """
    if hgccoord.obstime:
        obstime = hgccoord.obstime
    else:
        obstime = hgsframe.obstime
    s_lon = hgccoord.spherical.lon - _carrington_offset(obstime).to(
        u.deg)
    representation = SphericalWrap180Representation(s_lon, hgccoord.lat,
                                                    hgccoord.radius)
    return hgsframe.realize_frame(representation)


@frame_transform_graph.transform(FunctionTransform, Heliocentric,
                                 Helioprojective)
def hcc_to_hpc(helioccoord, heliopframe):
    """
    Convert from Heliocentic Cartesian to Helioprojective Cartesian.
    """
    # Propagate obstime explicitly.
    if heliopframe.obstime is None:
        heliopframe._obstime = helioccoord._obstime

    x = helioccoord.x.to(u.m)
    y = helioccoord.y.to(u.m)
    z = helioccoord.z.to(u.m)

    # d is calculated as the distance between the points
    # (x,y,z) and (0,0,D0).
    distance = np.sqrt(x**2 + y**2 + (helioccoord.observer.radius - z)**2)

    hpcx = np.rad2deg(np.arctan2(x, helioccoord.observer.radius - z))
    hpcy = np.rad2deg(np.arcsin(y / distance))

    representation = SphericalWrap180Representation(hpcx, hpcy,
                                                    distance.to(u.km))
    return heliopframe.realize_frame(representation)


@frame_transform_graph.transform(FunctionTransform, Helioprojective,
                                 Heliocentric)
def hpc_to_hcc(heliopcoord, heliocframe):
    """
    Convert from Helioprojective Cartesian to Heliocentric Cartesian.
    """
    if not isinstance(heliopcoord.observer, BaseCoordinateFrame):
        raise ConvertError("Cannot transform helioprojective coordinates to "
                           "heliocentric coordinates for observer '{}' "
                           "without `obstime` being specified.".format(heliopcoord.observer))

    heliopcoord = heliopcoord.calculate_distance()
    x = np.deg2rad(heliopcoord.Tx)
    y = np.deg2rad(heliopcoord.Ty)

    cosx = np.cos(x)
    sinx = np.sin(x)
    cosy = np.cos(y)
    siny = np.sin(y)

    rx = (heliopcoord.distance.to(u.m)) * cosy * sinx
    ry = (heliopcoord.distance.to(u.m)) * siny
    rz = (heliopcoord.observer.radius.to(u.m)) - (
        heliopcoord.distance.to(u.m)) * cosy * cosx

    representation = CartesianRepresentation(
        rx.to(u.km), ry.to(u.km), rz.to(u.km))
    return heliocframe.realize_frame(representation)


@frame_transform_graph.transform(FunctionTransform, Heliocentric,
                                 HeliographicStonyhurst)
def hcc_to_hgs(helioccoord, heliogframe):
    """
    Convert from Heliocentric Cartesian to Heliographic Stonyhurst.
    """
    if not isinstance(helioccoord.observer, BaseCoordinateFrame):
        raise ConvertError("Cannot transform heliocentric coordinates to "
                           "heliographic coordinates for observer '{}' "
                           "without `obstime` being specified.".format(helioccoord.observer))

    x = helioccoord.x.to(u.m)
    y = helioccoord.y.to(u.m)
    z = helioccoord.z.to(u.m)

    l0_rad = helioccoord.observer.lon
    b0_deg = helioccoord.observer.lat

    cosb = np.cos(np.deg2rad(b0_deg))
    sinb = np.sin(np.deg2rad(b0_deg))

    hecr = np.sqrt(x**2 + y**2 + z**2)
    hgln = np.arctan2(x, z * cosb - y * sinb) + l0_rad
    hglt = np.arcsin((y * cosb + z * sinb) / hecr)

    representation = SphericalWrap180Representation(
        np.rad2deg(hgln), np.rad2deg(hglt), hecr.to(u.km))
    return heliogframe.realize_frame(representation)


@frame_transform_graph.transform(FunctionTransform, HeliographicStonyhurst,
                                 Heliocentric)
def hgs_to_hcc(heliogcoord, heliocframe):
    """
    Convert from Heliographic Stonyhurst to Heliograpic Carrington.
    """
    hglon = heliogcoord.lon
    hglat = heliogcoord.lat
    r = heliogcoord.radius.to(u.m)

    if heliocframe.obstime is None:
        heliocframe._obstime = heliogcoord.obstime

    if not isinstance(heliocframe.observer, BaseCoordinateFrame):
        raise ConvertError("Cannot transform heliographic coordinates to "
                           "heliocentric coordinates for observer '{}' "
                           "without `obstime` being specified.".format(heliocframe.observer))

    l0_rad = heliocframe.observer.lon.to(u.rad)
    b0_deg = heliocframe.observer.lat

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

    representation = CartesianRepresentation(
        x.to(u.km), y.to(u.km), zz.to(u.km))

    return heliocframe.realize_frame(representation)


@frame_transform_graph.transform(FunctionTransform, Helioprojective,
                                 Helioprojective)
def hpc_to_hpc(heliopcoord, heliopframe):
    """
    This converts from HPC to HPC, with different observer location parameters.
    It does this by transforming through HGS.
    """
    if (heliopcoord.observer == heliopframe.observer or
        (heliopcoord.observer.lat == heliopframe.observer.lat and
         heliopcoord.observer.lon == heliopframe.observer.lon and
         heliopcoord.observer.radius == heliopframe.observer.radius)):
        return heliopframe.realize_frame(heliopcoord._data)

    if not isinstance(heliopframe.observer, BaseCoordinateFrame):
        raise ConvertError("Cannot transform between helioprojective frames "
                           "without `obstime` being specified for observer {}.".format(heliopframe.observer))
    if not isinstance(heliopcoord.observer, BaseCoordinateFrame):
        raise ConvertError("Cannot transform between helioprojective frames "
                           "without `obstime` being specified for observer {}.".format(heliopcoord.observer))

    hgs = heliopcoord.transform_to(HeliographicStonyhurst)
    hgs.observer = heliopframe.observer
    hpc = hgs.transform_to(heliopframe)

    return hpc


def _make_rotation_matrix(start_representation, end_representation):
    """
    Return the matrix for the direct rotation from one representation to a second representation.
    The representations need not be normalized first.
    """
    A = start_representation.to_cartesian()
    B = end_representation.to_cartesian()
    rotation_axis = A.cross(B)
    rotation_angle = -np.arccos(A.dot(B) / (A.norm() * B.norm()))  # negation is required

    # This line works around some input/output quirks of Astropy's rotation_matrix()
    matrix = np.array(rotation_matrix(rotation_angle, rotation_axis.xyz.value.tolist()))
    return matrix


# The Sun's north pole is oriented RA=286.13 deg, dec=63.87 deg in ICRS, and thus HCRS as well
# (See Archinal et al. 2011,
#   "Report of the IAU Working Group on Cartographic Coordinates and Rotational Elements: 2009")
# The orientation of the north pole in ICRS/HCRS is assumed to be constant in time
_solar_north_pole_hcrs = UnitSphericalRepresentation(lon=286.13*u.deg, lat=63.87*u.deg)


# Calculate the rotation matrix to de-tilt the Sun's rotation axis to be parallel to the Z axis
_sun_detilt_matrix = _make_rotation_matrix(_solar_north_pole_hcrs,
                                           CartesianRepresentation(0, 0, 1))


@frame_transform_graph.transform(DynamicMatrixTransform, HCRS, HeliographicStonyhurst)
def hcrs_to_hgs(hcrscoord, hgsframe):
    """
    Convert from HCRS to Heliographic Stonyhurst.
    """
    if hgsframe.obstime is None:
        raise ValueError("To perform this transformation the coordinate"
                         " Frame needs an obstime Attribute")

    # Determine the Sun-Earth vector in this de-tilt frame
    sun_pos_icrs = get_body_barycentric('sun', hgsframe.obstime)
    earth_pos_icrs = get_body_barycentric('earth', hgsframe.obstime)
    sun_earth_hcrs = earth_pos_icrs - sun_pos_icrs
    sun_earth_detilt = sun_earth_hcrs.transform(_sun_detilt_matrix)

    # Remove the component of the Sun-Earth vector that is parallel to the Sun's north pole
    hgs_x_axis_detilt = CartesianRepresentation(sun_earth_detilt.xyz * [1, 1, 0])

    # The above vector, which is in the Sun's equatorial plane, is also the X axis of HGS
    x_axis = CartesianRepresentation(1, 0, 0)
    rot_matrix = _make_rotation_matrix(hgs_x_axis_detilt, x_axis)

    return matrix_product(rot_matrix, _sun_detilt_matrix)


@frame_transform_graph.transform(DynamicMatrixTransform, HeliographicStonyhurst, HCRS)
def hgs_to_hcrs(hgscoord, hcrsframe):
    """
    Convert from Heliographic Stonyhurst to HCRS.
    """
    return matrix_transpose(hcrs_to_hgs(hcrsframe, hgscoord))


__doc__ += _make_transform_graph_docs()
