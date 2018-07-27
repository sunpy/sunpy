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
from astropy.coordinates import (HCRS, BaseCoordinateFrame, ConvertError,
                                 get_body_barycentric)
from astropy.coordinates.baseframe import frame_transform_graph
from astropy.coordinates.matrix_utilities import (matrix_product,
                                                  matrix_transpose,
                                                  rotation_matrix)
from astropy.coordinates.representation import (CartesianRepresentation,
                                                SphericalRepresentation,
                                                UnitSphericalRepresentation)
from astropy.coordinates.transformations import (DynamicMatrixTransform,
                                                 FunctionTransform)
from sunpy.sun import sun

from .frames import (Heliocentric, HeliographicCarrington,
                     HeliographicStonyhurst, Helioprojective)

try:
    from astropy.coordinates.builtin_frames import _make_transform_graph_docs as make_transform_graph_docs
except ImportError:
    from astropy.coordinates import make_transform_graph_docs as _make_transform_graph_docs
    make_transform_graph_docs = lambda: _make_transform_graph_docs(frame_transform_graph)

try:
    from astropy.units import allclose as quantity_allclose
except ImportError:
    from astropy.tests.helper import quantity_allclose


RSUN_METERS = sun.constants.get('radius').si.to(u.m)

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
    if hgcframe.obstime is None or np.any(hgcframe.obstime != hgscoord.obstime):
        raise ValueError("Can not transform from Heliographic Stonyhurst to "
                         "Heliographic Carrington, unless both frames have matching obstime.")

    c_lon = hgscoord.spherical.lon + _carrington_offset(hgscoord.obstime).to(u.deg)
    representation = SphericalRepresentation(c_lon, hgscoord.spherical.lat,
                                             hgscoord.spherical.distance)
    hgcframe = hgcframe.__class__(obstime=hgscoord.obstime)

    return hgcframe.realize_frame(representation)


@frame_transform_graph.transform(FunctionTransform, HeliographicCarrington,
                                 HeliographicStonyhurst)
def hgc_to_hgs(hgccoord, hgsframe):
    """
    Convert from Heliograpic Carrington to Heliographic Stonyhurst.
    """
    if hgsframe.obstime is None or np.any(hgsframe.obstime != hgccoord.obstime):
        raise ValueError("Can not transform from Heliographic Carrington to "
                         "Heliographic Stonyhurst, unless both frames have matching obstime.")
    obstime = hgsframe.obstime
    s_lon = hgccoord.spherical.lon - _carrington_offset(obstime).to(
        u.deg)
    representation = SphericalRepresentation(s_lon, hgccoord.spherical.lat,
                                             hgccoord.spherical.distance)

    return hgsframe.realize_frame(representation)


@frame_transform_graph.transform(FunctionTransform, Heliocentric,
                                 Helioprojective)
def hcc_to_hpc(helioccoord, heliopframe):
    """
    Convert from Heliocentic Cartesian to Helioprojective Cartesian.
    """

    x = helioccoord.x.to(u.m)
    y = helioccoord.y.to(u.m)
    z = helioccoord.z.to(u.m)

    # d is calculated as the distance between the points
    # (x,y,z) and (0,0,D0).
    distance = np.sqrt(x**2 + y**2 + (helioccoord.observer.radius - z)**2)

    hpcx = np.rad2deg(np.arctan2(x, helioccoord.observer.radius - z))
    hpcy = np.rad2deg(np.arcsin(y / distance))

    representation = SphericalRepresentation(hpcx, hpcy,
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

    representation = SphericalRepresentation(
        np.rad2deg(hgln), np.rad2deg(hglt), hecr.to(u.km))
    return heliogframe.realize_frame(representation)


@frame_transform_graph.transform(FunctionTransform, HeliographicStonyhurst,
                                 Heliocentric)
def hgs_to_hcc(heliogcoord, heliocframe):
    """
    Convert from Heliographic Stonyhurst to Heliocentric Cartesian.
    """
    hglon = heliogcoord.spherical.lon
    hglat = heliogcoord.spherical.lat
    r = heliogcoord.spherical.distance
    if r.unit is u.one and quantity_allclose(r, 1*u.one):
        r = np.ones_like(r)
        r *= RSUN_METERS

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
        (quantity_allclose(heliopcoord.observer.lat, heliopframe.observer.lat) and
         quantity_allclose(heliopcoord.observer.lon, heliopframe.observer.lon) and
         quantity_allclose(heliopcoord.observer.radius, heliopframe.observer.radius))):
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


def _make_rotation_matrix_from_reprs(start_representation, end_representation):
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
_SOLAR_NORTH_POLE_HCRS = UnitSphericalRepresentation(lon=286.13*u.deg, lat=63.87*u.deg)


# Calculate the rotation matrix to de-tilt the Sun's rotation axis to be parallel to the Z axis
_SUN_DETILT_MATRIX = _make_rotation_matrix_from_reprs(_SOLAR_NORTH_POLE_HCRS,
                                                      CartesianRepresentation(0, 0, 1))


@frame_transform_graph.transform(DynamicMatrixTransform, HCRS, HeliographicStonyhurst)
def hcrs_to_hgs(hcrscoord, hgsframe):
    """
    Convert from HCRS to Heliographic Stonyhurst (HGS).

    HGS shares the same origin (the Sun) as HCRS, but has its Z axis aligned with the Sun's
    rotation axis and its X axis aligned with the projection of the Sun-Earth vector onto the Sun's
    equatorial plane (i.e., the component of the Sun-Earth vector perpendicular to the Z axis).
    Thus, the transformation matrix is the product of the matrix to align the Z axis (by de-tilting
    the Sun's rotation axis) and the matrix to align the X axis.  The first matrix is independent
    of time and is pre-computed, while the second matrix depends on the time-varying Sun-Earth
    vector.
    """
    if hgsframe.obstime is None:
        raise ValueError("To perform this transformation the coordinate"
                         " Frame needs an obstime Attribute")

    # Determine the Sun-Earth vector in ICRS
    # Since HCRS is ICRS with an origin shift, this is also the Sun-Earth vector in HCRS
    sun_pos_icrs = get_body_barycentric('sun', hgsframe.obstime)
    earth_pos_icrs = get_body_barycentric('earth', hgsframe.obstime)
    sun_earth = earth_pos_icrs - sun_pos_icrs

    # De-tilt the Sun-Earth vector to the frame with the Sun's rotation axis parallel to the Z axis
    sun_earth_detilt = sun_earth.transform(_SUN_DETILT_MATRIX)

    # Remove the component of the Sun-Earth vector that is parallel to the Sun's north pole
    # (The additional transpose operations are to handle both scalar and array obstime situations)
    hgs_x_axis_detilt = CartesianRepresentation((sun_earth_detilt.xyz.T * [1, 1, 0]).T)

    # The above vector, which is in the Sun's equatorial plane, is also the X axis of HGS
    x_axis = CartesianRepresentation(1, 0, 0)
    if hgsframe.obstime.isscalar:
        rot_matrix = _make_rotation_matrix_from_reprs(hgs_x_axis_detilt, x_axis)
    else:
        rot_matrix_list = [_make_rotation_matrix_from_reprs(vect, x_axis) for vect in hgs_x_axis_detilt]
        rot_matrix = np.stack(rot_matrix_list)

    return matrix_product(rot_matrix, _SUN_DETILT_MATRIX)


@frame_transform_graph.transform(DynamicMatrixTransform, HeliographicStonyhurst, HCRS)
def hgs_to_hcrs(hgscoord, hcrsframe):
    """
    Convert from Heliographic Stonyhurst to HCRS.
    """
    return matrix_transpose(hcrs_to_hgs(hcrsframe, hgscoord))


@frame_transform_graph.transform(FunctionTransform, HeliographicStonyhurst, HeliographicStonyhurst)
def hgs_to_hgs(from_coo, to_frame):
    if np.all(from_coo.obstime == to_frame.obstime):
        return to_frame.realize_frame(from_coo.data)
    else:
        return from_coo.transform_to(HCRS).transform_to(to_frame)


__doc__ += make_transform_graph_docs()
