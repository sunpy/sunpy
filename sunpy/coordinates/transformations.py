# -*- coding: utf-8 -*-
"""
Coordinate Transformation Functions

This module contains the functions for converting one
`sunpy.coordinates.frames` object to another.
"""
from __future__ import absolute_import, division

import numpy as np

from astropy import units as u
from astropy.coordinates.representation import (CartesianRepresentation,
                                                UnitSphericalRepresentation)
from astropy.coordinates.baseframe import frame_transform_graph
from astropy.coordinates.transformations import FunctionTransform, DynamicMatrixTransform
from astropy.coordinates.matrix_utilities import rotation_matrix, matrix_product, matrix_transpose
from astropy.coordinates import HCRS, get_body_barycentric

from sunpy import sun

from .representation import (SphericalWrap180Representation,
                             UnitSphericalWrap180Representation)
from .frames import (HeliographicStonyhurst, HeliographicCarrington,
                     Heliocentric, Helioprojective)

__all__ = ['hgs_to_hgc', 'hgc_to_hgs', 'hcc_to_hpc',
           'hpc_to_hcc', 'hcc_to_hgs', 'hgs_to_hcc']


def _carrington_offset(dateobs):
    """
    Calculate the HG Longitude offest based on a time
    """
    if dateobs is None:
        raise ValueError("To perform this transformation the coordinate"
                         " Frame needs a dateobs Attribute")
    return sun.heliographic_solar_center(dateobs)[0]

# =============================================================================
# ------------------------- Transformation Framework --------------------------
# =============================================================================


@frame_transform_graph.transform(FunctionTransform, HeliographicStonyhurst,
                                 HeliographicCarrington)
def hgs_to_hgc(hgscoord, hgcframe):
    """
    Transform from Heliographic Stonyhurst to Heliograpic Carrington.
    """
    c_lon = hgscoord.spherical.lon + _carrington_offset(hgscoord.dateobs).to(
        u.deg)
    representation = SphericalWrap180Representation(c_lon, hgscoord.lat,
                                                    hgscoord.radius)
    hgcframe = hgcframe.__class__(dateobs=hgscoord.dateobs)
    return hgcframe.realize_frame(representation)


@frame_transform_graph.transform(FunctionTransform, HeliographicCarrington,
                                 HeliographicStonyhurst)
def hgc_to_hgs(hgccoord, hgsframe):
    """
    Convert from Heliograpic Carrington to Heliographic Stonyhurst.
    """
    if hgccoord.dateobs:
        dateobs = hgccoord.dateobs
    else:
        dateobs = hgsframe.dateobs
    s_lon = hgccoord.spherical.lon - _carrington_offset(dateobs).to(
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
    if (heliopcoord.observer.lat == heliopframe.observer.lat and
        heliopcoord.observer.lon == heliopframe.observer.lon and
        heliopcoord.observer.radius == heliopframe.observer.radius):
        return heliopframe.realize_frame(heliopcoord._data)

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
    rotation_axis = A.cross(B) / (A.norm() * B.norm())
    rotation_angle = -np.arcsin(rotation_axis.norm()) # negation is required

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
    # Use the time in the HGS frame unless it is not defined
    if hgsframe.dateobs is not None:
        time = hgsframe.dateobs
    else:
        time = hcrscoord.obstime

    # Determine the Sun-Earth vector in this de-tilt frame
    sun_pos_icrs = get_body_barycentric('sun', time)
    earth_pos_icrs = get_body_barycentric('earth', time)
    sun_earth_hcrs = earth_pos_icrs - sun_pos_icrs
    sun_earth_detilt = sun_earth_hcrs.transform(_sun_detilt_matrix)

    # Remove the component of the Sun-Earth vector that is parallel to the Sun's north pole
    z_axis = CartesianRepresentation(0, 0, 1)
    hgs_x_axis_detilt = sun_earth_detilt - z_axis * sun_earth_detilt.dot(z_axis)

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


# Make a transformation graph for the documentation, borrowed lovingly from
# Astropy.


def _make_transform_graph_docs():
    """
    Generates a string for use with the coordinate package's docstring
    to show the available transforms and coordinate systems
    """
    import inspect
    from textwrap import dedent
    from sunpy.extern import six
    from astropy.coordinates.baseframe import (BaseCoordinateFrame,
                                               frame_transform_graph)

    import copy
    f = copy.deepcopy(frame_transform_graph)
    for f1 in frame_transform_graph._graph.keys():
        if 'sunpy' not in str(f1):
            del f._graph[f1]
        else:
            for f2 in frame_transform_graph._graph[f1].keys():
                if 'sunpy' not in str(f2):
                    del f._graph[f1][f2]

    # TODO: Make this just show the SunPy Frames
    isclass = inspect.isclass
    coosys = [item
              for item in list(six.itervalues(globals()))
              if isclass(item) and issubclass(item, BaseCoordinateFrame)]
    graphstr = f.to_dot_graph(addnodes=coosys, priorities=False)

    docstr = """
    The diagram below shows all of the coordinate systems built into the
    `~astropy.coordinates` package, their aliases (useful for converting
    other coordinates to them using attribute-style access) and the
    pre-defined transformations between them.  The user is free to
    override any of these transformations by defining new transformations
    between these systems, but the pre-defined transformations should be
    sufficient for typical usage.

    .. graphviz::

    """

    return dedent(docstr) + '    ' + graphstr.replace('\n', '\n    ')


__doc__ += _make_transform_graph_docs()
