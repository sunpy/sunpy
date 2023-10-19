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

.. note::
    This module requires the optional dependency `~spiceypy.spiceypy` to be
    installed.

Notes
-----
* All transformations from one SPICE frame to another SPICE frame go through
  `~astropy.coordinates.ICRS` as the intermediate frame, even if the origin
  shift to/from the solar-system barycenter is unnatural.  This also means that
  it is not possible to transform a 2D coordinate between frames because there
  is always an origin shift.
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
from astropy.coordinates import ICRS, SkyCoord, frame_transform_graph
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


# Registry of the generated frame classes
spice_frame_classes = []


# Defined for future use
class _SpiceBaseCoordinateFrame(SunPyBaseCoordinateFrame):
    pass


def _convert_to_et(time):
    return (time - _ET_REF_EPOCH).to_value('s')


def _make_astropy_name(spice_frame_name):
    # Replace plus/minus characters in the SPICE frame name with lowercase 'p'/'n'
    return f"spice_{spice_frame_name.translate(str.maketrans('+-', 'pn'))}"


def _install_frame_by_id(frame_id):
    frame_name = spiceypy.frmnam(frame_id)
    astropy_frame_name = _make_astropy_name(frame_name)

    frame_center, class_num, _ = spiceypy.frinfo(frame_id)
    frame_center_name = spiceypy.bodc2n(frame_center)
    log.info(f"Installing {frame_name} {_CLASS_TYPES[class_num]} frame ({frame_id}) "
             f"as '{astropy_frame_name}'")

    spice_frame = type(astropy_frame_name, (_SpiceBaseCoordinateFrame,), {})
    # Force the capitalization pattern of lowercase "spice_" followed by uppercase SPICE frame name
    spice_frame.name = spice_frame.__name__

    @frame_transform_graph.transform(FunctionTransformWithFiniteDifference, ICRS, spice_frame)
    def icrs_to_spice(from_icrs_coord, to_spice_frame):
        et = _convert_to_et(to_spice_frame.obstime)
        # matrix needs to be contiguous (see https://github.com/astropy/astropy/issues/15503)
        matrix = np.ascontiguousarray(spiceypy.sxform('J2000', frame_name, et)[..., :3, :3])
        icrs_offset = spiceypy.spkpos(frame_center_name,
                                      et, 'J2000', 'NONE', 'SSB')[0] << u.km
        shifted_old_pos = from_icrs_coord.cartesian - CartesianRepresentation(icrs_offset.T)
        new_pos = shifted_old_pos.transform(matrix)
        return to_spice_frame.realize_frame(new_pos)

    @frame_transform_graph.transform(FunctionTransformWithFiniteDifference, spice_frame, ICRS)
    def spice_to_icrs(from_spice_coord, to_icrs_frame):
        et = _convert_to_et(from_spice_coord.obstime)
        # matrix needs to be contiguous (see https://github.com/astropy/astropy/issues/15503)
        matrix = np.ascontiguousarray(spiceypy.sxform(frame_name, 'J2000', et)[..., :3, :3])
        shifted_new_pos = from_spice_coord.cartesian.transform(matrix)
        icrs_offset = spiceypy.spkpos(frame_center_name,
                                      et, 'J2000', 'NONE', 'SSB')[0] << u.km
        new_pos = shifted_new_pos + CartesianRepresentation(icrs_offset.T)
        return to_icrs_frame.realize_frame(new_pos)

    frame_transform_graph._add_merged_transform(spice_frame, ICRS, spice_frame)

    spice_frame_classes.append(spice_frame)


def _uninstall_frame_by_class(spice_frame_class):
    frame_transform_graph.remove_transform(ICRS, spice_frame_class, None)
    frame_transform_graph.remove_transform(spice_frame_class, ICRS, None)
    frame_transform_graph.remove_transform(spice_frame_class, spice_frame_class, None)
    del spice_frame_class


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
    if spice_frame_classes:
        log.info(f"Removing {len(spice_frame_classes)} existing SPICE frame classes")
        for spice_frame_class in spice_frame_classes:
            _uninstall_frame_by_class(spice_frame_class)
        spice_frame_classes.clear()

    # Generate all SPICE frame classes
    for class_num in _CLASS_TYPES.keys():
        frames = spiceypy.kplfrm(class_num)
        for frame_id in frames:
            _install_frame_by_id(frame_id)


@add_common_docstring(**_variables_for_parse_time_docstring())
def get_body(body, time, *, spice_frame_name='J2000', observer=None):
    """
    Get the location of a body via SPICE.

    Parameters
    ----------
    body : `int`, `str`
        The NAIF body ID, or a string that is resolvable to a body ID
    time : {parse_time_types}
        Time to use in a parse_time-compatible format.
    spice_frame_name : `str`
        The SPICE frame to use for the returned coordinate.  Defaults to ``'J2000'``,
        which is equivalent to Astropy's `~astropy.coordinates.ICRS`.
    observer : `~astropy.coordinates.SkyCoord`
        If `None`, the returned coordinate is the instantaneous or “true” location.
        If not `None`, the returned coordinate is the astrometric location (i.e.,
        accounts for light travel time to the specified observer).
    """
    body_name = spiceypy.bodc2n(body) if isinstance(body, int) else body
    obstime = parse_time(time)
    et = _convert_to_et(obstime)

    frame_center = spiceypy.frinfo(spiceypy.namfrm(spice_frame_name))[0]

    if observer is None:
        pos = spiceypy.spkpos(body_name,
                              et,
                              spice_frame_name,
                              'NONE',
                              spiceypy.bodc2n(frame_center))[0] << u.km
    else:
        obspos = observer.icrs.cartesian.xyz.to_value('km')
        pos, lt = spiceypy.spkcpo(body_name,
                                  et,
                                  spice_frame_name,
                                  'OBSERVER',
                                  'CN',
                                  obspos,
                                  'SSB',
                                  'J2000')
        log.info(f"Apparent body location accounts for {lt:.2f} seconds of light travel time")

        pos = pos[:3] << u.km
        if spice_frame_name != 'J2000':
            shift = spiceypy.spkpos(spiceypy.bodc2n(frame_center),
                                    et,
                                    'J2000',
                                    'NONE',
                                    'SSB')[0]
            obspos -= shift

            matrix = spiceypy.pxform('J2000', spice_frame_name, _convert_to_et(obstime))
            obspos = matrix @ obspos
        pos += obspos << u.km

    frame_name = 'icrs' if spice_frame_name == 'J2000' else _make_astropy_name(spice_frame_name)

    return SkyCoord(CartesianRepresentation(pos.T), frame=frame_name, obstime=obstime)
