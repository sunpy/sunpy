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
import logging
from copy import deepcopy
from functools import wraps

import erfa
import numpy as np

import astropy.units as u
from astropy.constants import c as speed_of_light
from astropy.coordinates import (
    HCRS,
    ICRS,
    ITRS,
    BaseCoordinateFrame,
    ConvertError,
    HeliocentricMeanEcliptic,
    get_body_barycentric,
)
from astropy.coordinates.baseframe import frame_transform_graph
from astropy.coordinates.builtin_frames import make_transform_graph_docs
from astropy.coordinates.builtin_frames.utils import get_jd12
from astropy.coordinates.matrix_utilities import matrix_transpose, rotation_matrix
from astropy.coordinates.representation import (
    CartesianRepresentation,
    SphericalRepresentation,
    UnitSphericalRepresentation,
)
from astropy.coordinates.transformations import FunctionTransform, FunctionTransformWithFiniteDifference
from astropy.time import Time

from sunpy import log
from sunpy.sun import constants
from sunpy.util.decorators import sunpycontextmanager
from .frames import (
    _J2000,
    GeocentricEarthEquatorial,
    GeocentricSolarEcliptic,
    GeocentricSolarMagnetospheric,
    Geomagnetic,
    Heliocentric,
    HeliocentricEarthEcliptic,
    HeliocentricInertial,
    HeliographicCarrington,
    HeliographicStonyhurst,
    Helioprojective,
    HelioprojectiveRadial,
    SolarMagnetic,
)

RSUN_METERS = constants.get('radius').si.to(u.m)

__all__ = ['transform_with_sun_center',
           'propagate_with_solar_surface']


# Boolean flag for whether to ignore the motion of the center of the Sun in inertial space
_ignore_sun_motion = False


# If not None, the name of the differential-rotation model to use for any obstime change
_autoapply_diffrot = None


@sunpycontextmanager
def transform_with_sun_center():
    """
    Context manager for coordinate transformations to ignore the motion of the center of the Sun.

    Normally, coordinates refer to a point in inertial space (relative to the barycenter of the
    solar system). Transforming to a different observation time does not move the point at all,
    but rather only updates the coordinate representation as needed for the origin and axis
    orientations at the new observation time. However, the center of the Sun moves over time.
    Thus, for example, a coordinate that lies on the surface of the Sun at one observation time
    will not continue to lie on the surface of the Sun at other observation times.

    Under this context manager, transformations will instead move the coordinate over time to
    "follow" the translational motion of the center of Sun, thus maintaining the position of the
    coordinate relative to the center of the Sun.

    Notes
    -----
    This context manager accounts only for the motion of the center of the Sun, i.e.,
    translational motion. The motion of solar features due to any rotation of the Sun about its
    rotational axis is not accounted for.

    Due to the implementation approach, this context manager modifies transformations between only
    these five coordinate frames:
    `~sunpy.coordinates.frames.HeliographicStonyhurst`,
    `~sunpy.coordinates.frames.HeliographicCarrington`,
    `~sunpy.coordinates.frames.HeliocentricInertial`,
    `~sunpy.coordinates.frames.Heliocentric`, and
    `~sunpy.coordinates.frames.Helioprojective`.

    Examples
    --------
    >>> from astropy.coordinates import SkyCoord
    >>> from sunpy.coordinates import HeliographicStonyhurst, transform_with_sun_center
    >>> import astropy.units as u
    >>> start_frame = HeliographicStonyhurst(obstime="2001-01-01")
    >>> end_frame = HeliographicStonyhurst(obstime="2001-02-01")
    >>> sun_center = SkyCoord(0*u.deg, 0*u.deg, 0*u.AU, frame=start_frame)
    >>> sun_center
    <SkyCoord (HeliographicStonyhurst: obstime=2001-01-01T00:00:00.000, rsun=695700.0 km): (lon, lat, radius) in (deg, deg, AU)
        (0., 0., 0.)>
    >>> sun_center.transform_to(end_frame)  # transformations do not normally follow Sun center
    <SkyCoord (HeliographicStonyhurst: obstime=2001-02-01T00:00:00.000, rsun=695700.0 km): (lon, lat, radius) in (deg, deg, AU)
        (23.33174233, -5.96399877, 0.00027959)>
    >>> with transform_with_sun_center():
    ...     sun_center.transform_to(end_frame)  # now following Sun center
    <SkyCoord (HeliographicStonyhurst: obstime=2001-02-01T00:00:00.000, rsun=695700.0 km): (lon, lat, radius) in (deg, deg, AU)
        (0., 0., 0.)>
    """
    try:
        global _ignore_sun_motion

        old_ignore_sun_motion = _ignore_sun_motion  # nominally False

        if not old_ignore_sun_motion:
            log.debug("Ignoring the motion of the center of the Sun for transformations")
        _ignore_sun_motion = True
        yield
    finally:
        if not old_ignore_sun_motion:
            log.debug("Stop ignoring the motion of the center of the Sun for transformations")
        _ignore_sun_motion = old_ignore_sun_motion


@sunpycontextmanager
def propagate_with_solar_surface(rotation_model='howard'):
    """
    Context manager for coordinate transformations to automatically apply solar
    differential rotation for any change in observation time.

    Normally, coordinates refer to a point in inertial space (relative to the
    barycenter of the solar system). Transforming to a different observation time
    does not move the point at all, but rather only updates the coordinate
    representation as needed for the origin and axis orientations at the new
    observation time.

    Under this context manager, transformations will instead treat the coordinate
    as if it were referring to a point on the solar surface instead of a point in
    inertial space. If a transformation has a change in observation time, the
    heliographic longitude of the point will be updated according to the specified
    rotation model.

    Parameters
    ----------
    rotation_model : `str`
        Accepted model names are ``'howard'`` (default), ``'snodgrass'``,
        ``'allen'``, and ``'rigid'``. See the documentation for
        :func:`~sunpy.sun.models.differential_rotation` for the differences
        between these models.

    Notes
    -----
    This context manager also ignores the motion of the center of the Sun (see
    :func:`~sunpy.coordinates.transform_with_sun_center`).

    Due to the implementation approach, this context manager modifies
    transformations between only these five coordinate frames:
    `~sunpy.coordinates.frames.HeliographicStonyhurst`,
    `~sunpy.coordinates.frames.HeliographicCarrington`,
    `~sunpy.coordinates.frames.HeliocentricInertial`,
    `~sunpy.coordinates.frames.Heliocentric`, and
    `~sunpy.coordinates.frames.Helioprojective`.

    Examples
    --------
    .. minigallery:: sunpy.coordinates.propagate_with_solar_surface

    >>> import astropy.units as u
    >>> from astropy.coordinates import SkyCoord
    >>> from sunpy.coordinates import HeliocentricInertial, propagate_with_solar_surface
    >>> meridian = SkyCoord(0*u.deg, [-60, -30, 0, 30, 60]*u.deg, 1*u.AU,
    ...                     frame=HeliocentricInertial, obstime='2021-09-15')
    >>> out_frame = HeliocentricInertial(obstime='2021-09-21')
    >>> with propagate_with_solar_surface():
    ...     print(meridian.transform_to(out_frame))
    <SkyCoord (HeliocentricInertial: obstime=2021-09-21T00:00:00.000): (lon, lat, distance) in (deg, deg, AU)
        [(70.24182965, -60., 1.),
         (82.09298036, -30., 1.),
         (85.9579703 ,   0., 1.),
         (82.09298036,  30., 1.),
         (70.24182965,  60., 1.)]>
    >>> with propagate_with_solar_surface(rotation_model='rigid'):
    ...     print(meridian.transform_to(out_frame))
    <SkyCoord (HeliocentricInertial: obstime=2021-09-21T00:00:00.000): (lon, lat, distance) in (deg, deg, AU)
        [(85.1064, -60., 1.), (85.1064, -30., 1.),
         (85.1064,   0., 1.), (85.1064,  30., 1.),
         (85.1064,  60., 1.)]>
    """
    with transform_with_sun_center():
        try:
            global _autoapply_diffrot

            old_autoapply_diffrot = _autoapply_diffrot  # nominally False

            log.debug("Enabling automatic solar differential rotation "
                      f"('{rotation_model}') for any changes in obstime")
            _autoapply_diffrot = rotation_model
            yield
        finally:
            if not old_autoapply_diffrot:
                log.debug("Disabling automatic solar differential rotation "
                          "for any changes in obstime")
            _autoapply_diffrot = old_autoapply_diffrot


# Global counter to keep track of the layer of transformation
_layer_level = 0


def _transformation_debug(description):
    """
    Decorator to produce debugging output for a transformation function: its description, inputs,
    and output. Unicode box-drawing characters are used.
    """
    def decorator(func):
        @wraps(func)
        def wrapped_func(*args, **kwargs):
            global _layer_level

            # Check if the logging level is at least DEBUG (for performance reasons)
            debug_output = log.getEffectiveLevel() <= logging.DEBUG

            if debug_output:
                # Indentation for transformation layer
                indentation = "\u2502   " * _layer_level

                # For the input arguments, add indentation to any lines after the first line
                from_str = str(args[0]).replace("\n", f"\n       {indentation}\u2502       ")
                to_str = str(args[1]).replace("\n", f"\n       {indentation}\u2502       ")

                # Log the description and the input arguments
                log.debug(f"{indentation}{description}")
                log.debug(f"{indentation}\u251c\u2500From: {from_str}")
                log.debug(f"{indentation}\u251c\u2500To  : {to_str}")

                # Increment the layer level to increase the indentation for nested transformations
                _layer_level += 1

            result = func(*args, **kwargs)

            if debug_output:
                # Decrement the layer level
                _layer_level -= 1

                # For the output, add intention to any lines after the first line
                out_str = str(result).replace("\n", f"\n       {indentation}        ")

                # Log the output
                log.debug(f"{indentation}\u2514\u2500Out : {out_str}")

            return result
        return wrapped_func
    return decorator


def _observers_are_equal(obs_1, obs_2):
    # Note that this also lets pass the situation where both observers are None
    if obs_1 is obs_2:
        return True

    # obs_1 != obs_2
    if obs_1 is None:
        raise ConvertError("The source observer is set to None, but the transformation requires "
                           "the source observer to be specified, as the destination observer "
                           f"is set to {obs_2}.")
    if obs_2 is None:
        raise ConvertError("The destination observer is set to None, but the transformation "
                           "requires the destination observer to be specified, as the "
                           f"source observer is set to {obs_1}.")
    if isinstance(obs_1, str):
        if obs_1 == "self":
            return False
        raise ConvertError("The source observer needs to have `obstime` set because the "
                           "destination observer is different.")
    if isinstance(obs_2, str):
        if obs_2 == "self":
            return False
        raise ConvertError("The destination observer needs to have `obstime` set because the "
                           "source observer is different.")

    return np.atleast_1d(u.allclose(obs_1.lat, obs_2.lat) and
                         u.allclose(obs_1.lon, obs_2.lon) and
                         u.allclose(obs_1.radius, obs_2.radius) and
                         _times_are_equal(obs_1.obstime, obs_2.obstime)).all()


def _check_observer_defined(frame):
    if frame.observer is None:
        raise ConvertError("This transformation cannot be performed because the "
                           f"{frame.__class__.__name__} frame has observer=None.")
    elif isinstance(frame.observer, str):
        if frame.observer != "self":
            raise ConvertError("This transformation cannot be performed because the "
                               f"{frame.__class__.__name__} frame needs a specified obstime "
                               f"to fully resolve observer='{frame.observer}'.")
        elif not isinstance(frame, HeliographicCarrington):
            raise ConvertError(f"The {frame.__class__.__name__} frame has observer='self' "
                               "but this is valid for only HeliographicCarrington frames.")


def _times_are_equal(time_1, time_2):
    # Checks whether times are equal
    if isinstance(time_1, Time) and isinstance(time_2, Time):
        # We explicitly perform the check in TAI to avoid possible numerical precision differences
        # between a time in UTC and the same time after a UTC->TAI->UTC conversion
        return np.all(time_1.tai == time_2.tai)

    # We also deem the times equal if they are both None
    return time_1 is None and time_2 is None


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
    if obstime is None or _times_are_equal(frame.obstime, obstime):
        return frame

    # Transform to the new obstime using the appropriate loopback transformation
    new_frame = frame.replicate(obstime=obstime)
    if frame.obstime is not None:
        return frame.transform_to(new_frame)
    else:
        return new_frame


def _rotation_matrix_hgs_to_hgc(obstime, observer_distance_from_sun):
    """
    Return the rotation matrix from HGS to HGC at the same observation time
    """
    if obstime is None:
        raise ConvertError("To perform this transformation, the coordinate"
                           " frame needs a specified `obstime`.")

    # Import here to avoid a circular import
    from .sun import L0, earth_distance

    # Calculate the difference in light travel time if the observer is at a different distance from
    # the Sun than the Earth is
    delta_time = (observer_distance_from_sun - earth_distance(obstime)) / speed_of_light

    # Calculate the corresponding difference in apparent longitude
    delta_lon = delta_time * constants.sidereal_rotation_rate

    # Rotation is only in longitude, so only around the Z axis
    return rotation_matrix(-(L0(obstime) + delta_lon), 'z')


@frame_transform_graph.transform(FunctionTransformWithFiniteDifference,
                                 HeliographicStonyhurst, HeliographicCarrington)
@_transformation_debug("HGS->HGC")
def hgs_to_hgc(hgscoord, hgcframe):
    """
    Convert from Heliographic Stonyhurst to Heliographic Carrington.
    """
    _check_observer_defined(hgcframe)
    if isinstance(hgcframe.observer, str) and hgcframe.observer == "self":
        observer_radius = hgscoord.radius
    else:
        observer_radius = hgcframe.observer.radius

    # First transform the HGS coord to the HGC obstime
    int_coord = _transform_obstime(hgscoord, hgcframe.obstime)

    # Rotate from HGS to HGC
    total_matrix = _rotation_matrix_hgs_to_hgc(int_coord.obstime, observer_radius)
    newrepr = int_coord.cartesian.transform(total_matrix)

    return hgcframe._replicate(newrepr, obstime=int_coord.obstime)


@frame_transform_graph.transform(FunctionTransformWithFiniteDifference,
                                 HeliographicCarrington, HeliographicStonyhurst)
@_transformation_debug("HGC->HGS")
def hgc_to_hgs(hgccoord, hgsframe):
    """
    Convert from Heliographic Carrington to Heliographic Stonyhurst.
    """
    _check_observer_defined(hgccoord)

    hgccoord = hgccoord.make_3d()

    if isinstance(hgccoord.observer, str) and hgccoord.observer == "self":
        observer_radius = hgccoord.radius
    else:
        observer_radius = hgccoord.observer.radius

    # First transform the HGC coord to the HGS obstime
    int_coord = _transform_obstime(hgccoord, hgsframe.obstime)

    # Rotate from HGC to HGS
    total_matrix = matrix_transpose(_rotation_matrix_hgs_to_hgc(int_coord.obstime,
                                                                observer_radius))
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
@_transformation_debug("HCC->HPC")
def hcc_to_hpc(helioccoord, heliopframe):
    """
    Convert from Heliocentric Cartesian to Helioprojective Cartesian.
    """
    _check_observer_defined(helioccoord)
    _check_observer_defined(heliopframe)

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
@_transformation_debug("HPC->HCC")
def hpc_to_hcc(heliopcoord, heliocframe):
    """
    Convert from Helioprojective Cartesian to Heliocentric Cartesian.
    """
    _check_observer_defined(heliopcoord)
    _check_observer_defined(heliocframe)

    heliopcoord = heliopcoord.make_3d()

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
@_transformation_debug("HCC->HGS")
def hcc_to_hgs(helioccoord, heliogframe):
    """
    Convert from Heliocentric Cartesian to Heliographic Stonyhurst.
    """
    _check_observer_defined(helioccoord)

    # Transform the HCC observer (in HGS) to the HCC obstime in case it's different
    hcc_observer_at_hcc_obstime = _transform_obstime(helioccoord.observer, helioccoord.obstime)

    total_matrix = _rotation_matrix_hcc_to_hgs(hcc_observer_at_hcc_obstime.lon,
                                               hcc_observer_at_hcc_obstime.lat)

    # Transform from HCC to HGS at the HCC obstime
    newrepr = helioccoord.cartesian.transform(total_matrix)
    int_coord = HeliographicStonyhurst(newrepr, obstime=hcc_observer_at_hcc_obstime.obstime)

    # For historical reasons, we support HCC with no obstime transforming to HGS with an obstime
    if int_coord.obstime is None and heliogframe.obstime is not None:
        int_coord = int_coord.replicate(obstime=heliogframe.obstime)

    # Loopback transform HGS as needed
    return int_coord.transform_to(heliogframe)


@frame_transform_graph.transform(FunctionTransformWithFiniteDifference,
                                 HeliographicStonyhurst, Heliocentric)
@_transformation_debug("HGS->HCC")
def hgs_to_hcc(heliogcoord, heliocframe):
    """
    Convert from Heliographic Stonyhurst to Heliocentric Cartesian.
    """
    _check_observer_defined(heliocframe)

    heliogcoord = heliogcoord.make_3d()

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
@_transformation_debug("HPC->HPC")
def hpc_to_hpc(from_coo, to_frame):
    """
    This converts from HPC to HPC, with different observer location parameters.
    It does this by transforming through HGS.
    """
    if _observers_are_equal(from_coo.observer, to_frame.observer) and \
       _times_are_equal(from_coo.obstime, to_frame.obstime):
        return to_frame.realize_frame(from_coo.data)

    _check_observer_defined(from_coo)
    _check_observer_defined(to_frame)

    hgs = from_coo.transform_to(HeliographicStonyhurst(obstime=to_frame.obstime))
    hpc = hgs.transform_to(to_frame)

    return hpc


def _matrix_hpc_to_hpr():
    # Returns the transformation matrix that permutes/swaps axes from HPC to HPR

    # HPR spherical coordinates are a right-handed frame with these equivalent Cartesian axes:
    #   HPR_X = HPC_Z
    #   HPR_Y = -HPC_Y
    #   HPR_Z = -HPC_X
    # (HPC_X and HPC_Y are not to be confused with HPC_Tx and HPC_Ty)
    return np.array([[0, 0, 1],
                     [0, -1, 0],
                     [-1, 0, 0]])


@frame_transform_graph.transform(FunctionTransformWithFiniteDifference,
                                 Helioprojective, HelioprojectiveRadial)
@_transformation_debug("HPC->HPR")
def hpc_to_hpr(hpc_coord, hpr_frame):
    """
    Convert from Helioprojective Cartesian to Helioprojective Radial.
    """
    _check_observer_defined(hpc_coord)
    _check_observer_defined(hpr_frame)

    # Transform the HPR observer (in HGS) to the HPR obstime in case it's different
    observer = _transform_obstime(hpr_frame.observer, hpr_frame.obstime)

    # Loopback transform HPC coord to obstime and observer of HPR frame
    int_frame = Helioprojective(obstime=observer.obstime, observer=observer)
    int_coord = hpc_coord.transform_to(int_frame)

    # Permute/swap axes from HPC to HPR equivalent Cartesian
    newrepr = int_coord.cartesian.transform(_matrix_hpc_to_hpr())

    return hpr_frame.realize_frame(newrepr)


@frame_transform_graph.transform(FunctionTransformWithFiniteDifference,
                                 HelioprojectiveRadial, Helioprojective)
@_transformation_debug("HPR->HPC")
def hpr_to_hpc(hpr_coord, hpc_frame):
    """
    Convert from Helioprojective Radial to Helioprojective Cartesian.
    """
    _check_observer_defined(hpr_coord)
    _check_observer_defined(hpc_frame)

    # Permute/swap axes from HPR to HPC equivalent Cartesian
    newrepr = hpr_coord.cartesian.transform(matrix_transpose(_matrix_hpc_to_hpr()))

    # Transform the HPR observer (in HGS) to the HPR obstime in case it's different
    observer = _transform_obstime(hpr_coord.observer, hpr_coord.obstime)

    # Complete the conversion of HPR to HPC at the obstime and observer of the HPR coord
    int_coord = Helioprojective(newrepr, obstime=observer.obstime, observer=observer)

    # Loopback transform HPC as needed
    return int_coord.transform_to(hpc_frame)


frame_transform_graph._add_merged_transform(HelioprojectiveRadial, Helioprojective, HelioprojectiveRadial)


def _rotation_matrix_reprs_to_reprs(start_representation, end_representation):
    """
    Return the matrix for the direct rotation from one representation to a second representation.
    The representations need not be normalized first, and can be arrays of representations.
    """
    A = start_representation.to_cartesian()
    B = end_representation.to_cartesian()
    rotation_axis = A.cross(B)
    # Calculate the angle using both cross and dot products to minimize numerical-precision issues
    rotation_angle = -np.arctan2(rotation_axis.norm(), A.dot(B))  # negation is required

    if rotation_angle.isscalar:
        # This line works around some input/output quirks of Astropy's rotation_matrix()
        matrix = np.array(rotation_matrix(rotation_angle, rotation_axis.xyz.value.tolist()))
    else:
        matrix_list = [np.array(rotation_matrix(angle, axis.xyz.value.tolist()))
                       for angle, axis in zip(rotation_angle, rotation_axis)]
        matrix = np.stack(matrix_list)

    return matrix


def _rotation_matrix_reprs_to_xz_about_z(representations):
    """
    Return one or more matrices for rotating one or more representations around the Z axis into the
    XZ plane.
    """
    A = representations.to_cartesian()

    # Zero out the Z components
    # (The additional transpose operations are to handle both scalar and array inputs)
    A_no_z = CartesianRepresentation((A.xyz.T * [1, 1, 0]).T)

    # Rotate the resulting vector to the X axis
    x_axis = CartesianRepresentation(1, 0, 0)
    matrix = _rotation_matrix_reprs_to_reprs(A_no_z, x_axis)

    return matrix


def _sun_earth_icrf(time):
    """
    Return the Sun-Earth vector for ICRF-based frames.
    """
    sun_pos_icrs = get_body_barycentric('sun', time)
    earth_pos_icrs = get_body_barycentric('earth', time)
    return earth_pos_icrs - sun_pos_icrs


# The Sun's north pole is oriented RA=286.13 deg, dec=63.87 deg in ICRS, and thus HCRS as well
# (See Archinal et al. 2011,
#   "Report of the IAU Working Group on Cartographic Coordinates and Rotational Elements: 2009")
# The orientation of the north pole in ICRS/HCRS is assumed to be constant in time
_SOLAR_NORTH_POLE_HCRS = UnitSphericalRepresentation(lon=constants.get('alpha_0'),
                                                     lat=constants.get('delta_0'))


# Calculate the rotation matrix to de-tilt the Sun's rotation axis to be parallel to the Z axis
_SUN_DETILT_MATRIX = _rotation_matrix_reprs_to_reprs(_SOLAR_NORTH_POLE_HCRS,
                                                     CartesianRepresentation(0, 0, 1))


def _affine_params_hcrs_to_hgs(hcrs_time, hgs_time):
    """
    Return the affine parameters (matrix and offset) from HCRS to HGS

    HGS shares the same origin (the Sun) as HCRS, but has its Z axis aligned with the Sun's
    rotation axis and its X axis aligned with the projection of the Sun-Earth vector onto the Sun's
    equatorial plane (i.e., the component of the Sun-Earth vector perpendicular to the Z axis).
    Thus, the transformation matrix is the product of the matrix to align the Z axis (by de-tilting
    the Sun's rotation axis) and the matrix to align the X axis. The first matrix is independent
    of time and is pre-computed, while the second matrix depends on the time-varying Sun-Earth
    vector.
    """
    # Determine the Sun-Earth vector in ICRS
    # Since HCRS is ICRS with an origin shift, this is also the Sun-Earth vector in HCRS
    sun_pos_icrs = get_body_barycentric('sun', hgs_time)
    earth_pos_icrs = get_body_barycentric('earth', hgs_time)
    sun_earth = earth_pos_icrs - sun_pos_icrs

    # De-tilt the Sun-Earth vector to the frame with the Sun's rotation axis parallel to the Z axis
    sun_earth_detilt = sun_earth.transform(_SUN_DETILT_MATRIX)

    # Rotate the Sun-Earth vector about the Z axis so that it lies in the XZ plane
    rot_matrix = _rotation_matrix_reprs_to_xz_about_z(sun_earth_detilt)

    total_matrix = rot_matrix @ _SUN_DETILT_MATRIX

    # All of the above is calculated for the HGS observation time
    # If the HCRS observation time is different, calculate the translation in origin
    if not _ignore_sun_motion and np.any(hcrs_time != hgs_time):
        sun_pos_old_icrs = get_body_barycentric('sun', hcrs_time)
        offset_icrf = sun_pos_old_icrs - sun_pos_icrs
    else:
        offset_icrf = sun_pos_icrs * 0  # preserves obstime shape

    offset = offset_icrf.transform(total_matrix)
    return total_matrix, offset


@frame_transform_graph.transform(FunctionTransformWithFiniteDifference,
                                 HCRS, HeliographicStonyhurst)
@_transformation_debug("HCRS->HGS")
def hcrs_to_hgs(hcrscoord, hgsframe):
    """
    Convert from HCRS to Heliographic Stonyhurst (HGS).

    Even though we calculate the parameters for the affine transform, we use
    ``FunctionTransformWithFiniteDifference`` because otherwise there is no way to account for the
    induced angular velocity when transforming a coordinate with velocity information.
    """
    if hgsframe.obstime is None:
        raise ConvertError("To perform this transformation, the HeliographicStonyhurst"
                           " frame needs a specified `obstime`.")

    rot_matrix, offset = _affine_params_hcrs_to_hgs(hcrscoord.obstime, hgsframe.obstime)

    return hgsframe.realize_frame(hcrscoord.cartesian.transform(rot_matrix) + offset)


@frame_transform_graph.transform(FunctionTransformWithFiniteDifference,
                                 HeliographicStonyhurst, HCRS)
@_transformation_debug("HGS->HCRS")
def hgs_to_hcrs(hgscoord, hcrsframe):
    """
    Convert from Heliographic Stonyhurst to HCRS.

    Even though we calculate the parameters for the affine transform, we use
    ``FunctionTransformWithFiniteDifference`` because otherwise there is no way to account for the
    induced angular velocity when transforming a coordinate with velocity information.
    """
    if hgscoord.obstime is None:
        raise ConvertError("To perform this transformation, the HeliographicStonyhurst"
                           " frame needs a specified `obstime`.")

    hgscoord = hgscoord.make_3d()

    # Calculate the matrix and offset in the HCRS->HGS direction
    forward_matrix, forward_offset = _affine_params_hcrs_to_hgs(hcrsframe.obstime, hgscoord.obstime)

    # Invert the transformation to get the HGS->HCRS transformation
    reverse_matrix = matrix_transpose(forward_matrix)
    reverse_offset = (-forward_offset).transform(reverse_matrix)

    return hcrsframe.realize_frame(hgscoord.cartesian.transform(reverse_matrix) + reverse_offset)


@frame_transform_graph.transform(FunctionTransformWithFiniteDifference,
                                 HeliographicStonyhurst, HeliographicStonyhurst)
@_transformation_debug("HGS->HGS")
def hgs_to_hgs(from_coo, to_frame):
    """
    Convert between two Heliographic Stonyhurst frames.
    """
    if to_frame.obstime is None:
        return from_coo.replicate()
    elif _times_are_equal(from_coo.obstime, to_frame.obstime):
        return to_frame.realize_frame(from_coo.data)
    else:
        if _autoapply_diffrot:
            from_coo = from_coo._apply_diffrot((to_frame.obstime - from_coo.obstime).to('day'),
                                               _autoapply_diffrot)
        return from_coo.transform_to(HCRS(obstime=to_frame.obstime)).transform_to(to_frame)


@frame_transform_graph.transform(FunctionTransformWithFiniteDifference,
                                 HeliographicCarrington, HeliographicCarrington)
@_transformation_debug("HGC->HGC")
def hgc_to_hgc(from_coo, to_frame):
    """
    Convert between two Heliographic Carrington frames.
    """
    if _observers_are_equal(from_coo.observer, to_frame.observer) and \
       _times_are_equal(from_coo.obstime, to_frame.obstime):
        return to_frame.realize_frame(from_coo.data)

    _check_observer_defined(from_coo)
    _check_observer_defined(to_frame)

    # Convert through HGS
    hgscoord = from_coo.transform_to(HeliographicStonyhurst(obstime=from_coo.obstime))

    return hgscoord.transform_to(to_frame)


@frame_transform_graph.transform(FunctionTransformWithFiniteDifference,
                                 Heliocentric, Heliocentric)
@_transformation_debug("HCC->HCC")
def hcc_to_hcc(from_coo, to_frame):
    """
    Convert between two Heliocentric frames.
    """
    if _observers_are_equal(from_coo.observer, to_frame.observer) and \
       _times_are_equal(from_coo.obstime, to_frame.obstime):
        return to_frame.realize_frame(from_coo.data)

    _check_observer_defined(from_coo)
    _check_observer_defined(to_frame)

    # Convert through HGS
    hgscoord = from_coo.transform_to(HeliographicStonyhurst(obstime=to_frame.obstime))

    return hgscoord.transform_to(to_frame)


def _rotation_matrix_hme_to_hee(hmeframe):
    """
    Return the rotation matrix from HME to HEE at the same observation time
    """
    # Get the Sun-Earth vector
    sun_earth = HCRS(_sun_earth_icrf(hmeframe.obstime), obstime=hmeframe.obstime)
    sun_earth_hme = sun_earth.transform_to(hmeframe).cartesian

    # Rotate the Sun-Earth vector about the Z axis so that it lies in the XZ plane
    rot_matrix = _rotation_matrix_reprs_to_xz_about_z(sun_earth_hme)

    # Tilt the rotated Sun-Earth vector so that it is aligned with the X axis
    tilt_matrix = _rotation_matrix_reprs_to_reprs(sun_earth_hme.transform(rot_matrix),
                                                  CartesianRepresentation(1, 0, 0))

    return tilt_matrix @ rot_matrix


@frame_transform_graph.transform(FunctionTransformWithFiniteDifference,
                                 HeliocentricMeanEcliptic, HeliocentricEarthEcliptic)
@_transformation_debug("HME->HEE")
def hme_to_hee(hmecoord, heeframe):
    """
    Convert from Heliocentric Mean Ecliptic to Heliocentric Earth Ecliptic
    """
    if heeframe.obstime is None:
        raise ConvertError("To perform this transformation, the coordinate"
                           " frame needs a specified `obstime`.")

    # Convert to the HME frame with mean equinox of date at the HEE obstime, through HCRS
    int_frame = HeliocentricMeanEcliptic(obstime=heeframe.obstime, equinox=heeframe.obstime)
    int_coord = hmecoord.transform_to(HCRS(obstime=hmecoord.obstime)).transform_to(int_frame)

    # Rotate the intermediate coord to the HEE frame
    total_matrix = _rotation_matrix_hme_to_hee(int_frame)
    newrepr = int_coord.cartesian.transform(total_matrix)

    return heeframe.realize_frame(newrepr)


@frame_transform_graph.transform(FunctionTransformWithFiniteDifference,
                                 HeliocentricEarthEcliptic, HeliocentricMeanEcliptic)
@_transformation_debug("HEE->HME")
def hee_to_hme(heecoord, hmeframe):
    """
    Convert from Heliocentric Earth Ecliptic to Heliocentric Mean Ecliptic
    """
    if heecoord.obstime is None:
        raise ConvertError("To perform this transformation, the coordinate"
                           " frame needs a specified `obstime`.")

    int_frame = HeliocentricMeanEcliptic(obstime=heecoord.obstime, equinox=heecoord.obstime)

    # Rotate the HEE coord to the intermediate frame
    total_matrix = matrix_transpose(_rotation_matrix_hme_to_hee(int_frame))
    int_repr = heecoord.cartesian.transform(total_matrix)
    int_coord = int_frame.realize_frame(int_repr)

    # Convert to the HME frame through HCRS
    return int_coord.transform_to(HCRS(obstime=int_coord.obstime)).transform_to(hmeframe)


@frame_transform_graph.transform(FunctionTransformWithFiniteDifference,
                                 HeliocentricEarthEcliptic, HeliocentricEarthEcliptic)
@_transformation_debug("HEE->HEE")
def hee_to_hee(from_coo, to_frame):
    """
    Convert between two Heliocentric Earth Ecliptic frames.
    """
    if _times_are_equal(from_coo.obstime, to_frame.obstime):
        return to_frame.realize_frame(from_coo.data)
    elif to_frame.obstime is None:
        return from_coo
    else:
        return from_coo.transform_to(HCRS(obstime=from_coo.obstime)).transform_to(to_frame)


@frame_transform_graph.transform(FunctionTransformWithFiniteDifference,
                                 HeliocentricEarthEcliptic, GeocentricSolarEcliptic)
@_transformation_debug("HEE->GSE")
def hee_to_gse(heecoord, gseframe):
    """
    Convert from Heliocentric Earth Ecliptic to Geocentric Solar Ecliptic
    """
    # First transform the HEE coord to the GSE obstime
    int_coord = _transform_obstime(heecoord, gseframe.obstime)

    if int_coord.obstime is None:
        raise ConvertError("To perform this transformation, the coordinate"
                           " frame needs a specified `obstime`.")

    # Import here to avoid a circular import
    from .sun import earth_distance

    # Find the Earth-object vector in the intermediate frame
    sun_earth_int = earth_distance(int_coord.obstime) * CartesianRepresentation(1, 0, 0)
    earth_object_int = int_coord.cartesian - sun_earth_int

    # Flip the vector in X and Y, but leave Z untouched
    # (The additional transpose operations are to handle both scalar and array inputs)
    newrepr = CartesianRepresentation((earth_object_int.xyz.T * [-1, -1, 1]).T)

    return gseframe._replicate(newrepr, obstime=int_coord.obstime)


@frame_transform_graph.transform(FunctionTransformWithFiniteDifference,
                                 GeocentricSolarEcliptic, HeliocentricEarthEcliptic)
@_transformation_debug("GSE->HEE")
def gse_to_hee(gsecoord, heeframe):
    """
    Convert from Geocentric Solar Ecliptic to Heliocentric Earth Ecliptic
    """
    # First transform the GSE coord to the HEE obstime
    int_coord = _transform_obstime(gsecoord, heeframe.obstime)

    if int_coord.obstime is None:
        raise ConvertError("To perform this transformation, the coordinate"
                           " frame needs a specified `obstime`.")

    # Import here to avoid a circular import
    from .sun import earth_distance

    # Find the Sun-object vector in the intermediate frame
    earth_sun_int = earth_distance(int_coord.obstime) * CartesianRepresentation(1, 0, 0)
    sun_object_int = int_coord.cartesian - earth_sun_int

    # Flip the vector in X and Y, but leave Z untouched
    # (The additional transpose operations are to handle both scalar and array inputs)
    newrepr = CartesianRepresentation((sun_object_int.xyz.T * [-1, -1, 1]).T)

    return heeframe._replicate(newrepr, obstime=int_coord.obstime)


@frame_transform_graph.transform(FunctionTransformWithFiniteDifference,
                                 GeocentricSolarEcliptic, GeocentricSolarEcliptic)
@_transformation_debug("GSE->GSE")
def gse_to_gse(from_coo, to_frame):
    """
    Convert between two Geocentric Solar Ecliptic frames.
    """
    if _times_are_equal(from_coo.obstime, to_frame.obstime):
        return to_frame.realize_frame(from_coo.data)
    else:
        heecoord = from_coo.transform_to(HeliocentricEarthEcliptic(obstime=from_coo.obstime))
        return heecoord.transform_to(to_frame)


def _rotation_matrix_hgs_to_hci(obstime):
    """
    Return the rotation matrix from HGS to HCI at the same observation time
    """
    z_axis = CartesianRepresentation(0, 0, 1)*u.m
    if not obstime.isscalar:
        z_axis = z_axis._apply('repeat', obstime.size)

    # Get the ecliptic pole in HGS
    ecliptic_pole = HeliocentricMeanEcliptic(z_axis, obstime=obstime, equinox=_J2000)
    ecliptic_pole_hgs = ecliptic_pole.transform_to(HeliographicStonyhurst(obstime=obstime))

    # Rotate the ecliptic pole to the -YZ plane, which aligns the solar ascending node with the X
    # axis
    rot_matrix = _rotation_matrix_reprs_to_xz_about_z(ecliptic_pole_hgs.cartesian)
    xz_to_yz_matrix = rotation_matrix(-90*u.deg, 'z')

    return xz_to_yz_matrix @ rot_matrix


@frame_transform_graph.transform(FunctionTransformWithFiniteDifference,
                                 HeliographicStonyhurst, HeliocentricInertial)
@_transformation_debug("HGS->HCI")
def hgs_to_hci(hgscoord, hciframe):
    """
    Convert from Heliographic Stonyhurst to Heliocentric Inertial
    """

    hgscoord = hgscoord.make_3d()

    # First transform the HGS coord to the HCI obstime
    int_coord = _transform_obstime(hgscoord, hciframe.obstime)

    if int_coord.obstime is None:
        raise ConvertError("To perform this transformation, the coordinate"
                           " frame needs a specified `obstime`.")

    # Rotate from HGS to HCI
    total_matrix = _rotation_matrix_hgs_to_hci(int_coord.obstime)
    newrepr = int_coord.cartesian.transform(total_matrix)

    return hciframe._replicate(newrepr, obstime=int_coord.obstime)


@frame_transform_graph.transform(FunctionTransformWithFiniteDifference,
                                 HeliocentricInertial, HeliographicStonyhurst)
@_transformation_debug("HCI->HGS")
def hci_to_hgs(hcicoord, hgsframe):
    """
    Convert from Heliocentric Inertial to Heliographic Stonyhurst
    """
    # First transform the HCI coord to the HGS obstime
    int_coord = _transform_obstime(hcicoord, hgsframe.obstime)

    if int_coord.obstime is None:
        raise ConvertError("To perform this transformation, the coordinate"
                           " frame needs a specified `obstime`.")

    # Rotate from HCI to HGS
    total_matrix = matrix_transpose(_rotation_matrix_hgs_to_hci(int_coord.obstime))
    newrepr = int_coord.cartesian.transform(total_matrix)

    return hgsframe._replicate(newrepr, obstime=int_coord.obstime)


@frame_transform_graph.transform(FunctionTransformWithFiniteDifference,
                                 HeliocentricInertial, HeliocentricInertial)
@_transformation_debug("HCI->HCI")
def hci_to_hci(from_coo, to_frame):
    """
    Convert between two Heliocentric Inertial frames.
    """
    if _times_are_equal(from_coo.obstime, to_frame.obstime):
        return to_frame.realize_frame(from_coo.data)
    else:
        return from_coo.transform_to(HeliographicStonyhurst(obstime=from_coo.obstime)).\
            transform_to(to_frame)


def _rotation_matrix_obliquity(time):
    """
    Return the rotation matrix from Earth equatorial to ecliptic coordinates
    """
    return rotation_matrix(erfa.obl06(*get_jd12(time, 'tt'))*u.radian, 'x')


@frame_transform_graph.transform(FunctionTransformWithFiniteDifference,
                                 HeliocentricMeanEcliptic, GeocentricEarthEquatorial)
@_transformation_debug("HME->GEI")
def hme_to_gei(hmecoord, geiframe):
    """
    Convert from Heliocentric Mean Ecliptic to Geocentric Earth Equatorial
    """
    if geiframe.obstime is None:
        raise ConvertError("To perform this transformation, the coordinate"
                           " frame needs a specified `obstime`.")

    # Use an intermediate frame of HME at the GEI observation time, through HCRS
    int_frame = HeliocentricMeanEcliptic(obstime=geiframe.obstime, equinox=geiframe.equinox)
    int_coord = hmecoord.transform_to(HCRS(obstime=int_frame.obstime)).transform_to(int_frame)

    # Get the Sun-Earth vector in the intermediate frame
    sun_earth = HCRS(_sun_earth_icrf(int_frame.obstime), obstime=int_frame.obstime)
    sun_earth_int = sun_earth.transform_to(int_frame).cartesian

    # Find the Earth-object vector in the intermediate frame
    earth_object_int = int_coord.cartesian - sun_earth_int

    # Rotate from ecliptic to Earth equatorial
    rot_matrix = matrix_transpose(_rotation_matrix_obliquity(int_frame.equinox))
    newrepr = earth_object_int.transform(rot_matrix)

    return geiframe.realize_frame(newrepr)


@frame_transform_graph.transform(FunctionTransformWithFiniteDifference,
                                 GeocentricEarthEquatorial, HeliocentricMeanEcliptic)
@_transformation_debug("GEI->HME")
def gei_to_hme(geicoord, hmeframe):
    """
    Convert from Geocentric Earth Equatorial to Heliocentric Mean Ecliptic
    """
    if geicoord.obstime is None:
        raise ConvertError("To perform this transformation, the coordinate"
                           " frame needs a specified `obstime`.")

    # Use an intermediate frame of HME at the GEI observation time
    int_frame = HeliocentricMeanEcliptic(obstime=geicoord.obstime, equinox=geicoord.equinox)

    # Get the Sun-Earth vector in the intermediate frame
    sun_earth = HCRS(_sun_earth_icrf(int_frame.obstime), obstime=int_frame.obstime)
    sun_earth_int = sun_earth.transform_to(int_frame).cartesian

    # Rotate from Earth equatorial to ecliptic
    rot_matrix = _rotation_matrix_obliquity(int_frame.equinox)
    earth_object_int = geicoord.cartesian.transform(rot_matrix)

    # Find the Sun-object vector in the intermediate frame
    sun_object_int = earth_object_int + sun_earth_int  # add in this order to preserve the original units
    int_coord = int_frame.realize_frame(sun_object_int)

    # Convert to the final frame through HCRS
    return int_coord.transform_to(HCRS(obstime=int_coord.obstime)).transform_to(hmeframe)


@frame_transform_graph.transform(FunctionTransformWithFiniteDifference,
                                 GeocentricEarthEquatorial, GeocentricEarthEquatorial)
@_transformation_debug("GEI->GEI")
def gei_to_gei(from_coo, to_frame):
    """
    Convert between two Geocentric Earth Equatorial frames.
    """
    if _times_are_equal(from_coo.equinox, to_frame.equinox) and \
       _times_are_equal(from_coo.obstime, to_frame.obstime):
        return to_frame.realize_frame(from_coo.data)
    else:
        return from_coo.transform_to(HCRS(obstime=from_coo.obstime)).transform_to(to_frame)


@frame_transform_graph.transform(FunctionTransformWithFiniteDifference,
                                 ITRS, Geomagnetic)
@_transformation_debug("GEO->MAG")
def geo_to_mag(geocoord, magframe):
    """
    Convert from Geographic (GEO) to Geomagnetic (MAG)
    """
    if magframe.obstime is None:
        raise ConvertError("To perform this transformation, the coordinate"
                           " frame needs a specified `obstime`.")

    # First transform the GEO coord to the MAG obstime
    int_coord = _transform_obstime(geocoord, magframe.obstime)

    lon, lat = magframe.dipole_lonlat

    lat_matrix = rotation_matrix(90*u.deg - lat, 'y')
    lon_matrix = rotation_matrix(lon, 'z')
    rot_matrix = lat_matrix @ lon_matrix

    newrepr = int_coord.cartesian.transform(rot_matrix)

    return magframe.realize_frame(newrepr)


@frame_transform_graph.transform(FunctionTransformWithFiniteDifference,
                                 Geomagnetic, ITRS)
@_transformation_debug("MAG->GEO")
def mag_to_geo(magcoord, geoframe):
    """
    Convert from Geomagnetic (MAG) to Geographic (GEO)
    """
    if magcoord.obstime is None:
        raise ConvertError("To perform this transformation, the coordinate"
                           " frame needs a specified `obstime`.")

    lon, lat = magcoord.dipole_lonlat

    lat_matrix = rotation_matrix(90*u.deg - lat, 'y')
    lon_matrix = rotation_matrix(lon, 'z')
    rot_matrix = lat_matrix @ lon_matrix

    newrepr = magcoord.cartesian.transform(matrix_transpose(rot_matrix))
    int_frame = geoframe.replicate_without_data(obstime=magcoord.obstime)
    int_coord = int_frame.realize_frame(newrepr)

    return int_coord.transform_to(geoframe)


frame_transform_graph._add_merged_transform(Geomagnetic, ITRS, Geomagnetic)


@frame_transform_graph.transform(FunctionTransformWithFiniteDifference,
                                 Geomagnetic, SolarMagnetic)
@_transformation_debug("MAG->SM")
def mag_to_sm(magcoord, smframe):
    """
    Convert from Geomagnetic (MAG) to Solar Magnetic (SM)
    """
    if magcoord.obstime is None:
        raise ConvertError("To perform this transformation, the coordinate"
                           " frame needs a specified `obstime`.")

    # First transform the MAG coord to the SM obstime
    int_coord = _transform_obstime(magcoord, smframe.obstime)

    sun = HCRS(0*u.deg, 0*u.deg, 0*u.AU, obstime=int_coord.obstime).transform_to(int_coord)

    rot_matrix = _rotation_matrix_reprs_to_xz_about_z(sun.cartesian)

    newrepr = int_coord.cartesian.transform(rot_matrix)

    return smframe.realize_frame(newrepr)


@frame_transform_graph.transform(FunctionTransformWithFiniteDifference,
                                 SolarMagnetic, Geomagnetic)
@_transformation_debug("SM->MAG")
def sm_to_mag(smcoord, magframe):
    """
    Convert from Solar Magnetic (SM) to Geomagnetic (MAG)
    """
    if smcoord.obstime is None:
        raise ConvertError("To perform this transformation, the coordinate"
                           " frame needs a specified `obstime`.")

    int_frame = magframe.replicate_without_data(obstime=smcoord.obstime)

    sun = HCRS(0*u.deg, 0*u.deg, 0*u.AU, obstime=int_frame.obstime).transform_to(int_frame)

    rot_matrix = _rotation_matrix_reprs_to_xz_about_z(sun.cartesian)

    newrepr = smcoord.cartesian.transform(matrix_transpose(rot_matrix))
    int_coord = int_frame.realize_frame(newrepr)

    return int_coord.transform_to(magframe)


frame_transform_graph._add_merged_transform(SolarMagnetic, Geomagnetic, SolarMagnetic)


@frame_transform_graph.transform(FunctionTransformWithFiniteDifference,
                                 SolarMagnetic, GeocentricSolarMagnetospheric)
@_transformation_debug("SM->GSM")
def sm_to_gsm(smcoord, gsmframe):
    """
    Convert from Solar Magnetic (SM) to Geocentric Solar Magnetospheric (GSM)
    """
    if smcoord.obstime is None:
        raise ConvertError("To perform this transformation, the coordinate"
                           " frame needs a specified `obstime`.")

    # First transform the SM coord to the GSM obstime
    int_coord = _transform_obstime(smcoord, gsmframe.obstime)

    sun = HCRS(0*u.deg, 0*u.deg, 0*u.AU, obstime=int_coord.obstime).transform_to(int_coord)

    rot_matrix = rotation_matrix(-sun.spherical.lat, 'y')

    newrepr = int_coord.cartesian.transform(rot_matrix)

    return gsmframe.realize_frame(newrepr)


@frame_transform_graph.transform(FunctionTransformWithFiniteDifference,
                                 GeocentricSolarMagnetospheric, SolarMagnetic)
@_transformation_debug("GSM->SM")
def gsm_to_sm(gsmcoord, smframe):
    """
    Convert from Geocentric Solar Magnetospheric (GSM) to Solar Magnetic (SM)
    """
    if gsmcoord.obstime is None:
        raise ConvertError("To perform this transformation, the coordinate"
                           " frame needs a specified `obstime`.")

    int_frame = smframe.replicate_without_data(obstime=gsmcoord.obstime)

    sun = HCRS(0*u.deg, 0*u.deg, 0*u.AU, obstime=int_frame.obstime).transform_to(int_frame)

    rot_matrix = rotation_matrix(-sun.spherical.lat, 'y')

    newrepr = gsmcoord.cartesian.transform(matrix_transpose(rot_matrix))
    int_coord = int_frame.realize_frame(newrepr)

    return int_coord.transform_to(smframe)


frame_transform_graph._add_merged_transform(GeocentricSolarMagnetospheric, SolarMagnetic, GeocentricSolarMagnetospheric)


def _make_sunpy_graph():
    """
    Culls down the full transformation graph for SunPy purposes and returns the string version
    """
    # Frames to keep in the transformation graph
    keep_list = ['icrs', 'hcrs', 'heliocentrictrueecliptic', 'heliocentricmeanecliptic',
                 'heliographic_stonyhurst', 'heliographic_carrington',
                 'heliocentric', 'helioprojective', 'helioprojectiveradial',
                 'heliocentricearthecliptic', 'geocentricsolarecliptic',
                 'heliocentricinertial', 'geocentricearthequatorial',
                 'geomagnetic', 'solarmagnetic', 'geocentricsolarmagnetospheric',
                 'gcrs', 'precessedgeocentric', 'geocentrictrueecliptic', 'geocentricmeanecliptic',
                 'cirs', 'altaz', 'itrs']

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

    docstr = make_transform_graph_docs(small_graph)

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
                    'Heliocentric', 'Helioprojective', 'HelioprojectiveRadial',
                    'HeliocentricEarthEcliptic', 'GeocentricSolarEcliptic',
                    'HeliocentricInertial', 'GeocentricEarthEquatorial',
                    'Geomagnetic', 'SolarMagnetic', 'GeocentricSolarMagnetospheric']
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
