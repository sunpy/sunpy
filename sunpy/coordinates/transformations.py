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
from copy import deepcopy

import numpy as np

import astropy.units as u
from astropy.coordinates import (ICRS, HCRS, ConvertError, BaseCoordinateFrame,
                                 get_body_barycentric, get_body_barycentric_posvel)
from astropy.coordinates.baseframe import frame_transform_graph
from astropy.coordinates.representation import (CartesianRepresentation, SphericalRepresentation,
                                                UnitSphericalRepresentation, CartesianDifferential)
from astropy.coordinates.transformations import (FunctionTransform, AffineTransform,
                                                 FunctionTransformWithFiniteDifference)
from astropy.coordinates.matrix_utilities import rotation_matrix, matrix_transpose

from sunpy.sun import constants

from .frames import Heliocentric, Helioprojective, HeliographicCarrington, HeliographicStonyhurst

try:
    from astropy.coordinates.builtin_frames import _make_transform_graph_docs as make_transform_graph_docs
except ImportError:
    from astropy.coordinates import make_transform_graph_docs as _make_transform_graph_docs
    make_transform_graph_docs = lambda: _make_transform_graph_docs(frame_transform_graph)


RSUN_METERS = constants.get('radius').si.to(u.m)

__all__ = ['hgs_to_hgc', 'hgc_to_hgs', 'hcc_to_hpc',
           'hpc_to_hcc', 'hcc_to_hgs', 'hgs_to_hcc',
           'hpc_to_hpc',
           'hcrs_to_hgs', 'hgs_to_hcrs',
           'hgs_to_hgs', 'hgc_to_hgc', 'hcc_to_hcc']


def _observers_are_equal(obs_1, obs_2, string_ok=False):
    if string_ok:
        if obs_1 == obs_2:
            return True
    if not (isinstance(obs_1, BaseCoordinateFrame) and isinstance(obs_2, BaseCoordinateFrame)):
        raise ValueError("To compare two observers, both must be instances of BaseCoordinateFrame. "
                         "Cannot compare two observers {} and {}.".format(obs_1, obs_2))
    return np.atleast_1d((u.allclose(obs_1.lat, obs_2.lat) and
                          u.allclose(obs_1.lon, obs_2.lon) and
                          u.allclose(obs_1.radius, obs_2.radius) and
                          obs_1.obstime == obs_2.obstime)).all()


# =============================================================================
# ------------------------- Transformation Framework --------------------------
# =============================================================================


def _transform_obstime(frame, obstime):
    """
    Transform a frame to a new obstime using the appropriate loopback transformation.
    If the new obstime is None, no transformation is performed.
    If the frame's obstime is None, the frame is copied with the new obstime.
    """
    # If obstime is None or the obstime matches, nothing needs to be done
    if obstime is None or np.all(frame.obstime == obstime):
        return frame

    # Transform to the new obstime using the appropriate loopback transformation
    new_frame = frame.replicate(obstime=obstime)
    if frame.obstime is not None:
        return frame.transform_to(new_frame)
    else:
        return new_frame


def _rotation_matrix_hgs_to_hgc(obstime):
    """
    Return the rotation matrix from HGS to HGC at the same observation time
    """
    if obstime is None:
        raise ValueError("To perform this transformation the coordinate"
                         " Frame needs an obstime Attribute")

    # Import here to avoid a circular import
    from .sun import L0

    # Rotation is only in longitude, so only around the Z axis
    return rotation_matrix(-L0(obstime), 'z')


@frame_transform_graph.transform(FunctionTransformWithFiniteDifference,
                                 HeliographicStonyhurst, HeliographicCarrington)
def hgs_to_hgc(hgscoord, hgcframe):
    """
    Convert from Heliographic Stonyhurst to Heliographic Carrington.
    """
    # First transform the HGS coord to the HGC obstime
    int_coord = _transform_obstime(hgscoord, hgcframe.obstime)

    # Rotate from HGS to HGC
    total_matrix = _rotation_matrix_hgs_to_hgc(int_coord.obstime)
    newrepr = int_coord.cartesian.transform(total_matrix)

    return hgcframe._replicate(newrepr, obstime=int_coord.obstime)


@frame_transform_graph.transform(FunctionTransformWithFiniteDifference,
                                 HeliographicCarrington, HeliographicStonyhurst)
def hgc_to_hgs(hgccoord, hgsframe):
    """
    Convert from Heliographic Carrington to Heliographic Stonyhurst.
    """
    # First transform the HGC coord to the HGS obstime
    int_coord = _transform_obstime(hgccoord, hgsframe.obstime)

    # Rotate from HGC to HGS
    total_matrix = matrix_transpose(_rotation_matrix_hgs_to_hgc(int_coord.obstime))
    newrepr = int_coord.cartesian.transform(total_matrix)

    return hgsframe._replicate(newrepr, obstime=int_coord.obstime)


def _matrix_hcc_to_hpc():
    # Returns the transformation matrix that permutes/swaps axes from HCC to HPC

    # HPC spherical coordinates are a left-handed frame with these equivalent Cartesian axes:
    #   HPC_X = -HCC_Z
    #   HPC_Y = HCC_X
    #   HPC_Z = HCC_Y
    # (HPC_X and HPC_Y are not to be confused with HPC_Tx and HPC_Ty)
    return np.array([[0, 0, -1],
                     [1, 0, 0],
                     [0, 1, 0]])


@frame_transform_graph.transform(FunctionTransformWithFiniteDifference,
                                 Heliocentric, Helioprojective)
def hcc_to_hpc(helioccoord, heliopframe):
    """
    Convert from Heliocentric Cartesian to Helioprojective Cartesian.
    """
    # Transform the HPC observer (in HGS) to the HPC obstime in case it's different
    observer = _transform_obstime(heliopframe.observer, heliopframe.obstime)

    # Loopback transform HCC coord to obstime and observer of HPC frame
    int_frame = Heliocentric(obstime=observer.obstime, observer=observer)
    int_coord = helioccoord.transform_to(int_frame)

    # Shift the origin from the Sun to the observer
    distance = int_coord.observer.radius
    newrepr = int_coord.cartesian - CartesianRepresentation(0*u.m, 0*u.m, distance)

    # Permute/swap axes from HCC to HPC equivalent Cartesian
    newrepr = newrepr.transform(_matrix_hcc_to_hpc())

    # Explicitly represent as spherical because external code (e.g., wcsaxes) expects it
    return heliopframe.realize_frame(newrepr.represent_as(SphericalRepresentation))


@frame_transform_graph.transform(FunctionTransformWithFiniteDifference,
                                 Helioprojective, Heliocentric)
def hpc_to_hcc(heliopcoord, heliocframe):
    """
    Convert from Helioprojective Cartesian to Heliocentric Cartesian.
    """
    if not isinstance(heliopcoord.observer, BaseCoordinateFrame):
        raise ConvertError("Cannot transform helioprojective coordinates to "
                           "heliocentric coordinates for observer '{}' "
                           "without `obstime` being specified.".format(heliopcoord.observer))

    heliopcoord = heliopcoord.calculate_distance()

    # Permute/swap axes from HPC equivalent Cartesian to HCC
    newrepr = heliopcoord.cartesian.transform(matrix_transpose(_matrix_hcc_to_hpc()))

    # Transform the HPC observer (in HGS) to the HPC obstime in case it's different
    observer = _transform_obstime(heliopcoord.observer, heliopcoord.obstime)

    # Shift the origin from the observer to the Sun
    distance = observer.radius
    newrepr += CartesianRepresentation(0*u.m, 0*u.m, distance)

    # Complete the conversion of HPC to HCC at the obstime and observer of the HPC coord
    int_coord = Heliocentric(newrepr, obstime=observer.obstime, observer=observer)

    # Loopback transform HCC as needed
    return int_coord.transform_to(heliocframe)


def _rotation_matrix_hcc_to_hgs(longitude, latitude):
    # Returns the rotation matrix from HCC to HGS based on the observer longitude and latitude

    # Permute the axes of HCC to match HGS Cartesian equivalent
    #   HGS_X = HCC_Z
    #   HGS_Y = HCC_X
    #   HGS_Z = HCC_Y
    axes_matrix = np.array([[0, 0, 1],
                            [1, 0, 0],
                            [0, 1, 0]])

    # Rotate in latitude and longitude (sign difference because of direction difference)
    lat_matrix = rotation_matrix(latitude, 'y')
    lon_matrix = rotation_matrix(-longitude, 'z')

    return lon_matrix @ lat_matrix @ axes_matrix


@frame_transform_graph.transform(FunctionTransformWithFiniteDifference,
                                 Heliocentric, HeliographicStonyhurst)
def hcc_to_hgs(helioccoord, heliogframe):
    """
    Convert from Heliocentric Cartesian to Heliographic Stonyhurst.
    """
    if not isinstance(helioccoord.observer, BaseCoordinateFrame):
        raise ConvertError("Cannot transform heliocentric coordinates to "
                           "heliographic coordinates for observer '{}' "
                           "without `obstime` being specified.".format(helioccoord.observer))

    # Transform the HCC observer (in HGS) to the HCC obstime in case it's different
    hcc_observer_at_hcc_obstime = _transform_obstime(helioccoord.observer, helioccoord.obstime)

    total_matrix = _rotation_matrix_hcc_to_hgs(hcc_observer_at_hcc_obstime.lon,
                                               hcc_observer_at_hcc_obstime.lat)

    # Transform from HCC to HGS at the HCC obstime
    newrepr = helioccoord.cartesian.transform(total_matrix)
    int_coord = HeliographicStonyhurst(newrepr, obstime=hcc_observer_at_hcc_obstime.obstime)

    # Loopback transform HGS if there is a change in obstime
    return _transform_obstime(int_coord, heliogframe.obstime)


@frame_transform_graph.transform(FunctionTransformWithFiniteDifference,
                                 HeliographicStonyhurst, Heliocentric)
def hgs_to_hcc(heliogcoord, heliocframe):
    """
    Convert from Heliographic Stonyhurst to Heliocentric Cartesian.
    """
    if not isinstance(heliocframe.observer, BaseCoordinateFrame):
        raise ConvertError("Cannot transform heliographic coordinates to "
                           "heliocentric coordinates for observer '{}' "
                           "without `obstime` being specified.".format(heliocframe.observer))

    # Loopback transform HGS if there is a change in obstime
    int_coord = _transform_obstime(heliogcoord, heliocframe.obstime)

    # Transform the HCC observer (in HGS) to the HCC obstime in case it's different
    hcc_observer_at_hcc_obstime = _transform_obstime(heliocframe.observer, int_coord.obstime)

    total_matrix = matrix_transpose(_rotation_matrix_hcc_to_hgs(hcc_observer_at_hcc_obstime.lon,
                                                                hcc_observer_at_hcc_obstime.lat))

    # Transform from HGS to HCC at the same obstime
    newrepr = int_coord.cartesian.transform(total_matrix)
    return heliocframe.realize_frame(newrepr)


@frame_transform_graph.transform(FunctionTransformWithFiniteDifference,
                                 Helioprojective, Helioprojective)
def hpc_to_hpc(from_coo, to_frame):
    """
    This converts from HPC to HPC, with different observer location parameters.
    It does this by transforming through HGS.
    """
    if _observers_are_equal(from_coo.observer, to_frame.observer, string_ok=True) and \
       np.all(from_coo.obstime == to_frame.obstime):
        return to_frame.realize_frame(from_coo.data)

    if not isinstance(to_frame.observer, BaseCoordinateFrame):
        raise ConvertError("Cannot transform between helioprojective frames "
                           "without `obstime` being specified for observer {}.".format(to_frame.observer))
    if not isinstance(from_coo.observer, BaseCoordinateFrame):
        raise ConvertError("Cannot transform between helioprojective frames "
                           "without `obstime` being specified for observer {}.".format(from_coo.observer))

    hgs = from_coo.transform_to(HeliographicStonyhurst(obstime=to_frame.obstime))
    hpc = hgs.transform_to(to_frame)

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


@frame_transform_graph.transform(AffineTransform, HCRS, HeliographicStonyhurst)
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

    # Check whether differentials are involved on either end
    has_differentials = ((hcrscoord._data is not None and hcrscoord.data.differentials) or
                         (hgsframe._data is not None and hgsframe.data.differentials))

    # Determine the Sun-Earth vector in ICRS
    # Since HCRS is ICRS with an origin shift, this is also the Sun-Earth vector in HCRS
    # If differentials exist, also obtain Sun and Earth velocities
    if has_differentials:
        sun_pos_icrs, sun_vel = get_body_barycentric_posvel('sun', hgsframe.obstime)
        earth_pos_icrs, earth_vel = get_body_barycentric_posvel('earth', hgsframe.obstime)
    else:
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

    total_matrix = rot_matrix @ _SUN_DETILT_MATRIX

    # All of the above is calculated for the HGS observation time
    # If the HCRS observation time is different, calculate the translation in origin
    if np.any(hcrscoord.obstime != hgsframe.obstime):
        sun_pos_old_icrs = get_body_barycentric('sun', hcrscoord.obstime)
        offset_icrf = sun_pos_icrs - sun_pos_old_icrs
    else:
        offset_icrf = sun_pos_icrs * 0  # preserves obstime shape

    # Add velocity if needed (at the HGS observation time)
    if has_differentials:
        vel_icrf = (sun_vel - earth_vel).represent_as(CartesianDifferential)
        offset_icrf = offset_icrf.with_differentials(vel_icrf)

    offset = offset_icrf.transform(total_matrix)
    return total_matrix, offset


@frame_transform_graph.transform(AffineTransform, HeliographicStonyhurst, HCRS)
def hgs_to_hcrs(hgscoord, hcrsframe):
    """
    Convert from Heliographic Stonyhurst to HCRS.
    """
    # Calculate the matrix and offset in the HCRS->HGS direction
    total_matrix, offset = hcrs_to_hgs(hcrsframe, hgscoord)

    # Invert the transformation to get the HGS->HCRS transformation
    reverse_matrix = matrix_transpose(total_matrix)
    # If differentials exist, properly negate the velocity
    if offset.differentials:
        pos = -offset.without_differentials()
        vel = -offset.differentials['s']
        offset = pos.with_differentials(vel)
    else:
        offset = -offset
    reverse_offset = offset.transform(reverse_matrix)

    return reverse_matrix, reverse_offset


@frame_transform_graph.transform(FunctionTransformWithFiniteDifference,
                                 HeliographicStonyhurst, HeliographicStonyhurst)
def hgs_to_hgs(from_coo, to_frame):
    """
    Convert between two Heliographic Stonyhurst frames.
    """
    if to_frame.obstime is None:
        return from_coo.replicate()
    elif np.all(from_coo.obstime == to_frame.obstime):
        return to_frame.realize_frame(from_coo.data)
    else:
        return from_coo.transform_to(HCRS(obstime=from_coo.obstime)).transform_to(to_frame)


@frame_transform_graph.transform(FunctionTransformWithFiniteDifference,
                                 HeliographicCarrington, HeliographicCarrington)
def hgc_to_hgc(from_coo, to_frame):
    """
    Convert between two Heliographic Carrington frames.
    """
    if np.all(from_coo.obstime == to_frame.obstime):
        return to_frame.realize_frame(from_coo.data)
    else:
        return from_coo.transform_to(HeliographicStonyhurst(obstime=from_coo.obstime)).\
               transform_to(to_frame)


@frame_transform_graph.transform(FunctionTransformWithFiniteDifference,
                                 Heliocentric, Heliocentric)
def hcc_to_hcc(from_coo, to_frame):
    """
    Convert between two Heliocentric frames.
    """
    if _observers_are_equal(from_coo.observer, to_frame.observer, string_ok=True) and \
       np.all(from_coo.obstime == to_frame.obstime):
        return to_frame.realize_frame(from_coo.data)

    # Convert through HGS
    hgscoord = from_coo.transform_to(HeliographicStonyhurst(obstime=to_frame.obstime))

    return hgscoord.transform_to(to_frame)


def _make_sunpy_graph():
    """
    Culls down the full transformation graph for SunPy purposes and returns the string version
    """
    # Frames to keep in the transformation graph
    keep_list = ['icrs', 'hcrs', 'heliocentrictrueecliptic', 'heliocentricmeanecliptic',
                 'heliographic_stonyhurst', 'heliographic_carrington',
                 'heliocentric', 'helioprojective',
                 'gcrs', 'precessedgeocentric', 'geocentrictrueecliptic', 'geocentricmeanecliptic',
                 'cirs', 'altaz', 'itrs']

    global frame_transform_graph
    backup_graph = deepcopy(frame_transform_graph)

    small_graph = deepcopy(frame_transform_graph)
    cull_list = [name for name in small_graph.get_names() if name not in keep_list]
    cull_frames = [small_graph.lookup_name(name) for name in cull_list]

    for frame in cull_frames:
        # Remove the part of the graph where the unwanted frame is the source frame
        if frame in small_graph._graph:
            del small_graph._graph[frame]

        # Remove all instances of the unwanted frame as the destination frame
        for entry in small_graph._graph:
            if frame in small_graph._graph[entry]:
                del (small_graph._graph[entry])[frame]

    # Clean up the node list
    for name in cull_list:
        small_graph._cached_names.pop(name)

    _add_astropy_node(small_graph)

    # Overwrite the main transform graph
    frame_transform_graph = small_graph

    docstr = make_transform_graph_docs()

    # Restore the main transform graph
    frame_transform_graph = backup_graph

    # Make adjustments to the graph
    docstr = _tweak_graph(docstr)

    return docstr


def _add_astropy_node(graph):
    """
    Add an 'Astropy' node that links to an ICRS node in the graph
    """
    class Astropy(BaseCoordinateFrame):
        name = "REPLACE"

    @graph.transform(FunctionTransform, Astropy, ICRS)
    def fake_transform1():
        pass

    @graph.transform(FunctionTransform, ICRS, Astropy)
    def fake_transform2():
        pass


def _tweak_graph(docstr):
    # Remove Astropy's diagram description
    output = docstr[docstr.find('.. Wrap the graph'):]

    # Change the Astropy node
    output = output.replace('Astropy [shape=oval label="Astropy\\n`REPLACE`"]',
                            'Astropy [shape=box3d style=filled fillcolor=lightcyan '
                            'label="Other frames\\nin Astropy"]')

    # Change the Astropy<->ICRS links to black
    output = output.replace('ICRS -> Astropy[  color = "#783001" ]',
                            'ICRS -> Astropy[  color = "#000000" ]')
    output = output.replace('Astropy -> ICRS[  color = "#783001" ]',
                            'Astropy -> ICRS[  color = "#000000" ]')

    # Set the nodes to be filled and cyan by default
    output = output.replace('AstropyCoordinateTransformGraph {',
                            'AstropyCoordinateTransformGraph {\n'
                            '        node [style=filled fillcolor=lightcyan]')

    # Set the nodes for SunPy frames to be white
    sunpy_frames = ['HeliographicStonyhurst', 'HeliographicCarrington',
                    'Heliocentric', 'Helioprojective']
    for frame in sunpy_frames:
        output = output.replace(frame + ' [', frame + ' [fillcolor=white ')

    # Set the rank direction to be left->right (as opposed to top->bottom)
    # Force nodes for ICRS, HCRS, and "Other frames in Astropy" to be at the same rank
    output = output.replace('        overlap=false',
                            '        overlap=false\n'
                            '        rankdir=LR\n'
                            '        {rank=same; ICRS; HCRS; Astropy}')

    output = output.replace('<ul>\n\n',
                            '<ul>\n\n' +
                            _add_legend_row('SunPy frames', 'white') +
                            _add_legend_row('Astropy frames', 'lightcyan'))

    return output


def _add_legend_row(label, color):
    row = '        <li style="list-style: none;">\n'\
          '            <p style="font-size: 12px;line-height: 24px;font-weight: normal;'\
          'color: #848484;padding: 0;margin: 0;">\n'\
          '                <b>' + label + ':</b>\n'\
          '                    <span class="dot" style="height: 20px;width: 40px;'\
          'background-color: ' + color + ';border-radius: 50%;border: 1px solid black;'\
          'display: inline-block;"></span>\n'\
          '            </p>\n'\
          '        </li>\n\n\n'
    return row
