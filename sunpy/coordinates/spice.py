"""
Experimental module to use the SkyCoord API for SPICE computations

.. warning::
    This module is under development, so may be subject to significant change.

The `SPICE <https://naif.jpl.nasa.gov/naif/>`__ observation geometry information
system is being increasingly used by space missions to describe the locations of
spacecraft and the time-varying orientations of reference frames.
While the `~spiceypy.spiceypy` package provides a Python interface for
performing SPICE computations, its API is very different from that of
`~astropy.coordinates.SkyCoord`.

This module "wraps" `~spiceypy.spiceypy` functionality so that relevant SPICE
computations can be accessed using the `~astropy.coordinates.SkyCoord` API.
When loading a set of kernels, a frame class and corresponding transformations
are created for each SPICE frame.  One can also query the location of a body
as computed via SPICE.

See :ref:`sphx_glr_generated_gallery_units_and_coordinates_spice.py` for an
example of how to use this module.

.. note::
    This module requires the optional dependency `~spiceypy.spiceypy` to be
    installed.

Notes
-----
* 2D coordinates can be transformed only if the to/from frames have the same
  SPICE body ID as their centers.
* Transformations of velocities are not yet supported.
* SPICE frame names are rendered as uppercase, except for plus/minus characters,
  which are replaced with lowercase ``'p'``/``'n'`` characters.
"""

import numpy as np

try:
    import spiceypy
except ImportError:
    raise ImportError("This module requires the optional dependency `spiceypy`.")

import astropy.units as u
from astropy.coordinates import ICRS, ConvertError, SkyCoord, frame_transform_graph
from astropy.coordinates.representation import CartesianRepresentation
from astropy.coordinates.transformations import FunctionTransformWithFiniteDifference
from astropy.time import Time

from sunpy import log
from sunpy.coordinates import SunPyBaseCoordinateFrame
from sunpy.time import parse_time
from sunpy.time.time import _variables_for_parse_time_docstring
from sunpy.util.decorators import add_common_docstring

__all__ = ['get_body', 'initialize']


# Note that this epoch is very slightly different from the typical definition of J2000.0 (in TT)
_ET_REF_EPOCH = Time('J2000', scale='tdb')

_CLASS_TYPES = {2: 'PCK', 3: 'CK', 4: 'TK', 5: 'dynamic', 6: 'switch'}


# Registry of the generated frame classes and center classes
_frame_registry = {}
_center_registry = {'SOLAR SYSTEM BARYCENTER': ICRS}


# Defined for future use
class _SpiceBaseCoordinateFrame(SunPyBaseCoordinateFrame):
    pass


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
    center_cls = type(astropy_center_name, (_SpiceBaseCoordinateFrame,), {})

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

    frame_cls = type(astropy_frame_name, (_SpiceBaseCoordinateFrame,), {})
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
    order to load multiple sets of kernels.  However, there may be unexpected
    behavior if this function is called after the frame classes start being used.
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
        The SPICE frame name to use for the returned coordinate.  Defaults to
        ``'J2000'``, which is equivalent to Astropy's `~astropy.coordinates.ICRS`.
    observer : `~astropy.coordinates.SkyCoord`
        If `None`, the returned coordinate is the instantaneous or “true” location.
        If not `None`, the returned coordinate is the astrometric location (i.e.,
        accounts for light travel time to the specified observer).
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
