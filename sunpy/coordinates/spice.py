"""
Bridge module to use the SkyCoord API for SPICE computations.

.. note::
    This module requires the optional dependency `~spiceypy.spiceypy` to be
    installed.

The `SPICE <https://naif.jpl.nasa.gov/naif/>`__ observation geometry information
system is being increasingly used by space missions to describe the locations of
spacecraft and the time-varying orientations of reference frames.
While the `~spiceypy.spiceypy` package provides a Python interface for
performing SPICE computations, its API is very different from that of
`~astropy.coordinates.SkyCoord`.

This module "wraps" `~spiceypy.spiceypy` functionality so that relevant SPICE
computations can be accessed using the `~astropy.coordinates.SkyCoord` API.
When loading a set of kernels, a frame class and corresponding transformations
are created for each SPICE frame. One can also query the location of a body
as computed via SPICE or retrieve the field of view (FOV) of an instrument.

To facilitate the conversion of a SPICE-based coordinate to the built-in frames
in `sunpy.coordinates`, every SPICE-based coordinate has the method
:meth:`~sunpy.coordinates.spice.SpiceBaseCoordinateFrame.to_helioprojective()`.
This method returns a coordinate in the `~sunpy.coordinates.Helioprojective`
frame with the ``observer`` frame attribute appropriately set.

Be aware that when converting a SPICE-based coordinate to/from a built-in frame,
there can be small inconsistencies due to differing planetary ephemerides and
models for various orientations.

Notes
-----
* 2D coordinates can be transformed only if the to/from frames have the same
  SPICE body ID as their centers.
* Transformations of velocities are not yet supported.
* SPICE frame names are rendered as uppercase, except for plus/minus characters,
  which are replaced with lowercase ``'p'``/``'n'`` characters.

Examples
--------
.. minigallery:: sunpy.coordinates.spice.initialize
"""
# Developer notes:
# * We create a public SkyCoord frame for each SPICE frame that is defined in
#   the kernels, but this does not include built-in SPICE frames (e.g., inertial
#   frames or IAU_* frames). The user needs to manually install each built-in
#   SPICE frame that they want to use because otherwise there are simply too
#   many.
# * We also create a private SkyCoord frame for each unique SPICE frame center.
#   Each SPICE frame defines its center, and typically many frames share the
#   same center. By creating these private frames for frame centers, we can
#   transform 2D coordinates between frames that share the same center because
#   the origin does not change.
# * Any transformation that involves a change in frame center (including even
#   a change in the body ID that still maps to the same location) will be
#   treated as a change in origin, and the transformation is routed through
#   ICRS. ICRS is a safe frame to use because the SPICE built-in inertial
#   frame 'J2000' is ICRS, despite its name.

import numpy as np

try:
    import spiceypy
except ImportError:
    raise ImportError("This module requires the optional dependency `spiceypy`.")

import astropy.units as u
from astropy.coordinates import ICRS, ConvertError, SkyCoord, frame_transform_graph
from astropy.coordinates.matrix_utilities import rotation_matrix
from astropy.coordinates.representation import CartesianRepresentation, SphericalRepresentation
from astropy.coordinates.transformations import FunctionTransformWithFiniteDifference
from astropy.time import Time

from sunpy import log
from sunpy.coordinates import SunPyBaseCoordinateFrame
from sunpy.time import parse_time
from sunpy.time.time import _variables_for_parse_time_docstring
from sunpy.util.decorators import add_common_docstring

__all__ = ['SpiceBaseCoordinateFrame', 'get_body', 'get_fov', 'initialize', 'install_frame', 'get_rotation_matrix']


# Note that this epoch is very slightly different from the typical definition of J2000.0 (in TT)
_ET_REF_EPOCH = Time('J2000', scale='tdb')

_CLASS_TYPES = {1: 'inertial', 2: 'PCK', 3: 'CK', 4: 'TK', 5: 'dynamic', 6: 'switch'}


# Registry of the generated frame classes and center classes
_frame_registry = {}
_center_registry = {'SOLAR SYSTEM BARYCENTER': ICRS}


@add_common_docstring(**_variables_for_parse_time_docstring())
class SpiceBaseCoordinateFrame(SunPyBaseCoordinateFrame):
    """
    Base class for all frames generated to represent SPICE frames.

    This class is not intended to be used directly and has no transformations
    defined.

    Parameters
    ----------
    obstime : {parse_time_types}
        The time of the observation. This is used to determine the
        position of solar-system bodies (e.g., the Sun and the Earth) as
        needed to define the origin and orientation of the frame.
    """
    def __init_subclass__(cls, **kwargs):
        cls._frame_name = kwargs.pop('frame_name', None)
        cls._center_name = kwargs.pop('center_name', None)

        super().__init_subclass__(**kwargs)

        cls.__doc__ = (f"Coordinate frame for the SPICE frame '{cls._frame_name}'\n\n"
                       f"Origin: '{cls._center_name}'\n\n"
                       "Parameters\n----------\n"
                       f"obstime : {_variables_for_parse_time_docstring()['parse_time_types']}\n"
                       "    The time of the observation. This is used to determine the\n"
                       "    position of solar-system bodies (e.g., the Sun and the Earth) as\n"
                       "    needed to define the origin and orientation of the frame.\n")

    def to_helioprojective(self):
        """
        Convert this coordinate to the Helioprojective frame.

        The center of the frame center is used as the ``observer`` of the
        `~sunpy.coordinates.Helioprojective` frame.

        Examples
        --------
        .. minigallery:: sunpy.coordinates.spice.SpiceBaseCoordinate.to_helioprojective
        """
        et = _convert_to_et(self.obstime)

        # Get the matrix to rotate from the SPICE frame to heliographic coordinates
        # matrix needs to be contiguous (see https://github.com/astropy/astropy/issues/15503)
        frame_to_iau = np.ascontiguousarray(spiceypy.sxform(self._frame_name, 'IAU_SUN', et)[..., :3, :3])

        # Get the observer location in heliographic coordinates
        obs_iau = spiceypy.spkpos(self._center_name, et, 'IAU_SUN', 'NONE', 'SUN')[0] << u.km
        obs_iau = CartesianRepresentation(obs_iau.T).represent_as(SphericalRepresentation)

        # Construct the matrix to rotate from heliographic coordinates to HPC-like coordinates
        iau_to_hpc = rotation_matrix(-obs_iau.lat, axis='y') @ rotation_matrix(obs_iau.lon, axis='z')

        # To get to actual HPC coordinates, need to flip the X axis
        flip_x = np.array([[-1, 0, 0], [0, 1, 0], [0, 0, 1]])

        # Transform the data by all of the sequence of these matrices to get the HPC vector
        hpc_repr = self.cartesian.transform(flip_x @ iau_to_hpc @ frame_to_iau)

        # Get the observer location in ICRS, with obstime define to be able to transform to HGS
        obs_icrs = spiceypy.spkpos(self._center_name, et, 'J2000', 'NONE', 'SSB')[0] << u.km
        obs_sc = SkyCoord(CartesianRepresentation(obs_icrs.T), frame='icrs', obstime=self.obstime)

        # Construct the HPC coordinate from the vector and the observer
        out_sc = SkyCoord(hpc_repr, frame='helioprojective', obstime=self.obstime, observer=obs_sc)
        if _is_2d(hpc_repr):
            out_sc.representation_type = 'unitspherical'
        return out_sc


def _convert_to_et(time):
    return (time - _ET_REF_EPOCH).to_value('s')


def _astropy_frame_name(spice_frame_name):
    # Replace plus/minus characters in the SPICE frame name with lowercase 'p'/'n'
    return f"spice_{spice_frame_name.translate(str.maketrans('+-', 'pn'))}"


def _astropy_center_name(spice_center_id):
    # Use the center ID directly, changing a negative sign to 'n'
    return f"_spice_center_{str(spice_center_id).replace('-', 'n')}"


def _is_2d(data):
    return data.norm().unit is u.one and u.allclose(data.norm(), 1*u.one)


def _install_center_by_id(center_id):
    center_name = spiceypy.bodc2n(center_id)
    if center_name in _center_registry.keys():
        return

    log.info(f"Creating ICRF frame with {center_name} ({center_id}) origin")
    astropy_center_name = _astropy_center_name(center_id)
    center_cls = type(astropy_center_name, (SpiceBaseCoordinateFrame,), {},
                      frame_name=None, center_name=center_name)

    @frame_transform_graph.transform(FunctionTransformWithFiniteDifference, ICRS, center_cls)
    def icrs_to_shifted(from_icrs_coord, to_shifted_frame):
        if _is_2d(from_icrs_coord.data):
            raise ConvertError("Cannot transform a 2D coordinate due to a shift in origin.")
        icrs_offset = spiceypy.spkpos(center_name,
                                      _convert_to_et(to_shifted_frame.obstime),
                                      'J2000', 'NONE', 'SSB')[0] << u.km
        shifted_pos = from_icrs_coord.cartesian - CartesianRepresentation(icrs_offset.T)
        return to_shifted_frame.realize_frame(shifted_pos)

    @frame_transform_graph.transform(FunctionTransformWithFiniteDifference, center_cls, ICRS)
    def shifted_to_icrs(from_shifted_coord, to_icrs_frame):
        if _is_2d(from_shifted_coord.data):
            raise ConvertError("Cannot transform a 2D coordinate due to a shift in origin.")
        icrs_offset = spiceypy.spkpos(center_name,
                                      _convert_to_et(from_shifted_coord.obstime),
                                      'J2000', 'NONE', 'SSB')[0] << u.km
        icrs_pos = from_shifted_coord.cartesian + CartesianRepresentation(icrs_offset.T)
        return to_icrs_frame.realize_frame(icrs_pos)

    frame_transform_graph._add_merged_transform(center_cls, ICRS, center_cls)

    _center_registry[center_name] = center_cls


def _install_frame_by_id(frame_id):
    frame_name = spiceypy.frmnam(frame_id)
    astropy_frame_name = _astropy_frame_name(frame_name)

    center_id, class_num, _ = spiceypy.frinfo(frame_id)
    center_name = spiceypy.bodc2n(center_id)
    log.info(f"Installing {frame_name} {_CLASS_TYPES[class_num]} frame ({frame_id}) "
             f"as '{astropy_frame_name}'")

    # Create a center class (if needed)
    _install_center_by_id(center_id)
    center_cls = _center_registry[center_name]

    frame_cls = type(astropy_frame_name, (SpiceBaseCoordinateFrame,), {},
                     frame_name=frame_name, center_name=center_name)
    # Force the capitalization pattern of lowercase "spice_" followed by uppercase SPICE frame name
    frame_cls.name = frame_cls.__name__

    @frame_transform_graph.transform(FunctionTransformWithFiniteDifference, center_cls, frame_cls)
    def rotate_from_icrf(from_shifted_coord, to_spice_frame):
        et = _convert_to_et(to_spice_frame.obstime)
        # matrix needs to be contiguous (see https://github.com/astropy/astropy/issues/15503)
        matrix = np.ascontiguousarray(spiceypy.sxform('J2000', frame_name, et)[..., :3, :3])
        new_pos = from_shifted_coord.data.transform(matrix)
        return to_spice_frame.realize_frame(new_pos)

    @frame_transform_graph.transform(FunctionTransformWithFiniteDifference, frame_cls, center_cls)
    def rotate_to_icrf(from_spice_coord, to_shifted_frame):
        et = _convert_to_et(from_spice_coord.obstime)
        # matrix needs to be contiguous (see https://github.com/astropy/astropy/issues/15503)
        matrix = np.ascontiguousarray(spiceypy.sxform(frame_name, 'J2000', et)[..., :3, :3])
        shifted_pos = from_spice_coord.data.transform(matrix)
        return to_shifted_frame.realize_frame(shifted_pos)

    frame_transform_graph._add_merged_transform(frame_cls, center_cls, frame_cls)

    _frame_registry[frame_name] = (frame_cls, center_cls)


def _uninstall_frame_by_class(target_class, from_class):
    frame_transform_graph.remove_transform(target_class, target_class, None)
    frame_transform_graph.remove_transform(from_class, target_class, None)
    frame_transform_graph.remove_transform(target_class, from_class, None)
    del target_class


def initialize(kernels):
    """
    Load one more more SPICE kernels and create corresponding frame classes.

    Parameters
    ----------
    kernels : `str`, `list` of `str`
         One or more SPICE kernel files

    Notes
    -----
    If a kernel file is a meta-kernel, make sure that the relative paths therein
    are correct for the current working directory, which may not be the same as the
    location of the meta-kernel file.

    This function has simple support for being called multiple times in a session in
    order to load multiple sets of kernels. However, there may be unexpected
    behavior if this function is called after the frame classes start being used.

    Examples
    --------
    .. minigallery:: sunpy.coordinates.spice.initialize
    """
    if not isinstance(kernels, list):
        kernels = [kernels]
    # furnsh() needs path strings
    spiceypy.furnsh([str(kernel) for kernel in kernels])

    # Remove all existing SPICE frame classes
    global _frame_registry
    if _frame_registry:
        log.info(f"Removing {len(_frame_registry)} existing SPICE frame classes")
        for spice_frame_name in _frame_registry.keys():
            frame_cls, center_cls = _frame_registry[spice_frame_name]
            _uninstall_frame_by_class(frame_cls, center_cls)
        _frame_registry.clear()

    # Remove all non-default SPICE center classes
    global _center_registry
    if len(_center_registry) > 1:
        log.info(f"Removing {len(_center_registry) - 1} existing SPICE origin classes")
        for spice_center_name, center_cls in _center_registry.items():
            if center_cls != ICRS:
                _uninstall_frame_by_class(center_cls, ICRS)
        _center_registry = {'SOLAR SYSTEM BARYCENTER': ICRS}

    # Generate all SPICE frame classes
    for class_num in _CLASS_TYPES.keys():
        frames = spiceypy.kplfrm(class_num)
        for frame_id in frames:
            install_frame(frame_id)


def install_frame(spice_frame):
    """
    Install a specified SPICE frame.

    Installing a SPICE frame creates a corresponding frame class. All frames
    defined in the kernel pool are already automatically installed in the call
    to :func:`~sunpy.coordinates.spice.initialize`, so this function is used to
    manually install built-in frames, namely inertial or body-fixed (PCK) frames.
    Some common built-in frames include 'IAU_SUN', 'IAU_EARTH', and 'ITRF93'.

    Parameters
    ----------
    spice_frame : `str`, `int`
        The SPICE frame name or frame ID to be installed.

    Examples
    --------
    .. minigallery:: sunpy.coordinates.spice.install_frame
    """
    if isinstance(spice_frame, str):
        frame_name = spice_frame.upper()
        frame_id = spiceypy.namfrm(frame_name)
        if frame_id == 0:
            raise ValueError(f"{frame_name} is not a valid SPICE frame name.")
    else:
        frame_id = spice_frame
        frame_name = spiceypy.frmnam(frame_id)
        if frame_name == '':
            raise ValueError(f"{frame_id} is not a valid SPICE frame ID.")

    if frame_name in _frame_registry:
        log.info(f"{frame_name} (frame_id) has already been installed.")
    else:
        _install_frame_by_id(frame_id)


@add_common_docstring(**_variables_for_parse_time_docstring())
def get_body(body, time, *, spice_frame='J2000', observer=None):
    """
    Get the location of a body via SPICE.

    Parameters
    ----------
    body : `int`, `str`
        The NAIF body ID, or a string that is resolvable to a body ID
    time : {parse_time_types}
        Time to use in a parse_time-compatible format.
    spice_frame : `str`
        The SPICE frame name to use for the returned coordinate. Defaults to
        ``'J2000'``, which is equivalent to Astropy's `~astropy.coordinates.ICRS`.
    observer : `~astropy.coordinates.SkyCoord`
        If `None`, the returned coordinate is the instantaneous or “true” location.
        If not `None`, the returned coordinate is the astrometric location (i.e.,
        accounts for light travel time to the specified observer).

    Examples
    --------
    .. minigallery:: sunpy.coordinates.spice.get_body
    """
    body_name = spiceypy.bodc2n(body) if isinstance(body, int) else body
    obstime = parse_time(time)
    et = _convert_to_et(obstime)

    frame_center = spiceypy.frinfo(spiceypy.namfrm(spice_frame))[0]

    if observer is None:
        pos = spiceypy.spkpos(body_name,
                              et,
                              spice_frame,
                              'NONE',
                              spiceypy.bodc2n(frame_center))[0] << u.km
    else:
        obspos = observer.icrs.cartesian.xyz.to_value('km')
        pos, lt = spiceypy.spkcpo(body_name,
                                  et,
                                  spice_frame,
                                  'OBSERVER',
                                  'CN',
                                  obspos,
                                  'SSB',
                                  'J2000')
        log.info(f"Apparent body location accounts for {lt:.2f} seconds of light travel time")

        pos = pos[:3] << u.km
        if spice_frame != 'J2000':
            shift = spiceypy.spkpos(spiceypy.bodc2n(frame_center),
                                    et,
                                    'J2000',
                                    'NONE',
                                    'SSB')[0]
            obspos -= shift

            matrix = spiceypy.pxform('J2000', spice_frame, _convert_to_et(obstime))
            obspos = matrix @ obspos
        pos += obspos << u.km

    frame_name = 'icrs' if spice_frame == 'J2000' else _astropy_frame_name(spice_frame)

    return SkyCoord(CartesianRepresentation(pos.T), frame=frame_name, obstime=obstime)


@add_common_docstring(**_variables_for_parse_time_docstring())
def get_fov(instrument, time, *, resolution=100):
    """
    Get the field of view (FOV) for an instrument via SPICE.

    Rectangular and polygonal FOVs are represented by their vertices. Circular FOVs
    are approximated by a series of points. This function does not yet support
    elliptical FOVs.

    .. note::
        The FOV determined from SPICE kernels may not be as accurate as the FOV
        obtained from other sources of information, particularly if the instrument
        is an imager.

    Parameters
    ----------
    instrument : `int`, `str`
        The NAIF ID for the instrument, or a string that is resolvable to an
        instrument ID
    time : {parse_time_types}
        Time to use in a parse_time-compatible format.
    resolution : `int`
        Number of points to use for a circular FOV. Defaults to 100.

    Examples
    --------
    .. minigallery:: sunpy.coordinates.spice.get_fov
    """
    instrument_name = spiceypy.bodc2n(instrument) if isinstance(instrument, int) else instrument
    obstime = parse_time(time)

    fov_shape, spice_frame, boresight, num_vectors, vectors = spiceypy.getfvn(instrument_name, 1000)

    if fov_shape == "ELLIPSE":
        raise ValueError("Elliptical FOVs are not yet supported.")

    # If circular FOV, we have only one vector, so we have to rotate for the rest of the points
    if fov_shape == "CIRCLE":
        angles = np.arange(0, 1, 1/resolution) * 360*u.deg
        matrix = rotation_matrix(angles, axis=boresight)
        vectors = matrix @ vectors[0]
        num_vectors = resolution

    # The FOV vectors need to be replicated if obstime is a time array
    if not obstime.isscalar:
        vectors = np.broadcast_to(vectors[:, np.newaxis, :], (num_vectors, *obstime.shape, 3))
        obstime = Time(np.tile(obstime, num_vectors)).reshape((num_vectors, *obstime.shape)).T

    fov = SkyCoord(CartesianRepresentation(vectors.T),
                   frame=_frame_registry[spice_frame][0],
                   obstime=obstime,
                   representation_type='unitspherical')
    return fov


@add_common_docstring(**_variables_for_parse_time_docstring())
def get_rotation_matrix(source_frame, target_frame, from_time, to_time=None):
    """
    Get the rotation matrix between the orientations of two SPICE frames.

    Parameters
    ----------
    source_frame : `str`
        The source frame specified by SPICE frame name.

    target_frame : `str`
        The target frame specified by SPICE frame name.

    from_time : {parse_time_types}
        The time of the source frame.

    to_time : {parse_time_types}
        The time of the target frame. Defaults to ``None``, which means
        ``from_time`` is used.

    Returns
    -------
    `~numpy.ndarray`
        A 3x3 rotation matrix for the change in orientation.

    Examples
    --------
    >>> from sunpy.coordinates.spice import get_rotation_matrix
    >>> source_frame = "J2000"
    >>> target_frame = "Galactic"
    >>> from_time = '2001-01-01T00:00:00'
    >>> rotation_matrix = get_rotation_matrix(source_frame, target_frame, from_time)
    >>> rotation_matrix
    array([[-0.05487554, -0.8734371 , -0.48383499],
           [ 0.49410945, -0.44482959,  0.74698225],
           [-0.86766614, -0.19807639,  0.45598379]])

    This rotation matrix can be used to re-orient a vector (field), e.g.:

    >>> vec_components = [1, 0, 0] * u.T
    >>> transformed_matrix = rotation_matrix @ vec_components
    >>> transformed_matrix
    <Quantity [-0.05487554, 0.49410945, -0.86766614] T>
    """
    source_frame_spice = source_frame.upper()
    target_frame_spice = target_frame.upper()

    from_time = parse_time(from_time)
    to_time = parse_time(to_time) if to_time else from_time

    from_time_et = _convert_to_et(from_time)
    to_time_et = _convert_to_et(to_time)

    # First transformation: from source frame at from_time to J2000
    from_source_to_j2000 = spiceypy.sxform(source_frame_spice, "J2000", from_time_et)

    # Second transformation: from J2000 at to_time to target frame
    from_j2000_to_target = spiceypy.sxform("J2000", target_frame_spice, to_time_et)

    # Combine: source -> J2000 -> target
    combined_transform = spiceypy.mxmg(from_j2000_to_target, from_source_to_j2000)

    # Extract the rotation matrix (upper left 3x3 block)
    rotation_matrix = combined_transform[:3, :3]

    return rotation_matrix
