import numpy as np

import astropy.units as u
import astropy.wcs.utils
from astropy.coordinates import (
    ITRS,
    BaseCoordinateFrame,
    CartesianRepresentation,
    SkyCoord,
    SphericalRepresentation,
)
from astropy.wcs import WCS

from sunpy import log
from .frames import (
    BaseCoordinateFrame,
    Heliocentric,
    HeliographicCarrington,
    HeliographicStonyhurst,
    Helioprojective,
    SunPyBaseCoordinateFrame,
)

__all__ = ['solar_wcs_frame_mapping', 'solar_frame_to_wcs_mapping']

try:
    # TODO: Remove vendored version after Astropy 5.0
    from astropy.wcs.utils import obsgeo_to_frame
except ImportError:
    def obsgeo_to_frame(obsgeo, obstime):
        """
        Convert a WCS obsgeo property into an `~builtin_frames.ITRS` coordinate frame.

        Parameters
        ----------
        obsgeo : array-like
            A shape ``(6, )`` array representing ``OBSGEO-[XYZ], OBSGEO-[BLH]`` as
            returned by ``WCS.wcs.obsgeo``.

        obstime : time-like
            The time assiociated with the coordinate, will be passed to
            `~.builtin_frames.ITRS` as the obstime keyword.

        Returns
        -------
        `~.builtin_frames.ITRS`
            An `~.builtin_frames.ITRS` coordinate frame
            representing the coordinates.

        Notes
        -----

        The obsgeo array as accessed on a `.WCS` object is a length 6 numpy array
        where the first three elements are the coordinate in a cartesian
        representation and the second 3 are the coordinate in a spherical
        representation.

        This function priorities reading the cartesian coordinates, and will only
        read the spherical coordinates if the cartesian coordinates are either all
        zero or any of the cartesian coordinates are non-finite.

        In the case where both the spherical and cartesian coordinates have some
        non-finite values the spherical coordinates will be returned with the
        non-finite values included.

        """
        if (obsgeo is None
            or len(obsgeo) != 6
            or np.all(np.array(obsgeo) == 0)
            or np.all(~np.isfinite(obsgeo))
        ):  # NOQA
            raise ValueError(f"Can not parse the 'obsgeo' location ({obsgeo}). "
                             "obsgeo should be a length 6 non-zero, finite numpy array")

        # If the cartesian coords are zero or have NaNs in them use the spherical ones
        if np.all(obsgeo[:3] == 0) or np.any(~np.isfinite(obsgeo[:3])):
            data = SphericalRepresentation(*(obsgeo[3:] * (u.deg, u.deg, u.m)))

        # Otherwise we assume the cartesian ones are valid
        else:
            data = CartesianRepresentation(*obsgeo[:3] * u.m)

        return ITRS(data, obstime=obstime)


def solar_wcs_frame_mapping(wcs):
    """
    This function registers the coordinates frames to their FITS-WCS coordinate
    type values in the `astropy.wcs.utils.wcs_to_celestial_frame` registry.

    Parameters
    ----------
    wcs : astropy.wcs.WCS

    Returns
    -------
    astropy.coordinates.BaseCoordinateFrame
    """

    if hasattr(wcs, "coordinate_frame"):
        return wcs.coordinate_frame

    dateobs = wcs.wcs.dateobs or None

    # Get observer coordinate from the WCS auxillary information
    required_attrs = {HeliographicStonyhurst: ['hgln_obs', 'hglt_obs', 'dsun_obs'],
                      HeliographicCarrington: ['crln_obs', 'hglt_obs', 'dsun_obs']}

    # Get rsun from the WCS auxillary information
    rsun = wcs.wcs.aux.rsun_ref
    if rsun is not None:
        rsun *= u.m

    # TODO: remove these errors in sunpy 4.1
    bad_attrs = [f'.{attr}' for attr in ['rsun', 'heliographic_observer']
                 if hasattr(wcs, attr)]
    if len(bad_attrs):
        raise ValueError(f"The {' and '.join(bad_attrs)} attribute(s) on a WCS "
                         "are no longer supported.")

    observer = None
    for frame, attr_names in required_attrs.items():
        attrs = [getattr(wcs.wcs.aux, attr_name) for attr_name in attr_names]
        if all([attr is not None for attr in attrs]):
            kwargs = {'obstime': dateobs}
            if rsun is not None:
                kwargs['rsun'] = rsun
            if issubclass(frame, HeliographicCarrington):
                kwargs['observer'] = 'self'

            observer = frame(attrs[0] * u.deg,
                             attrs[1] * u.deg,
                             attrs[2] * u.m,
                             **kwargs)

    # Read the observer out of obsgeo for ground based observers
    if observer is None:
        try:
            observer = obsgeo_to_frame(wcs.wcs.obsgeo, dateobs)
            observer = SkyCoord(observer, rsun=rsun)
        except ValueError as e:
            # The helper function assumes you know the obsgeo coords you are
            # parsing are good, we are not sure, so catch the error.

            # This approach could lead to an invalid observer (i.e. one of the
            # coords being NaN), but only if the WCS has been constructed like that.
            log.debug(f"Could not parse obsgeo coordinates from WCS:\n{e}")

    # Collect all of the possible frame attributes, although some may be removed later
    frame_args = {'obstime': dateobs}
    if observer is not None:
        frame_args['observer'] = observer
    if rsun is not None:
        frame_args['rsun'] = rsun

    frame_class = _sunpy_frame_class_from_ctypes(wcs.wcs.ctype)

    if frame_class:
        if frame_class == HeliographicStonyhurst:
            frame_args.pop('observer', None)
        if frame_class == Heliocentric:
            frame_args.pop('rsun', None)

        return frame_class(**frame_args)


def _sunpy_frame_class_from_ctypes(ctypes):
    # Truncate the ctype to the first four letters
    ctypes = {c[:4] for c in ctypes}

    mapping = {
        Helioprojective: {'HPLN', 'HPLT'},
        HeliographicStonyhurst: {'HGLN', 'HGLT'},
        HeliographicCarrington: {'CRLN', 'CRLT'},
        Heliocentric: {'SOLX', 'SOLY'},
    }

    for frame_class, ctype_pair in mapping.items():
        if ctype_pair <= ctypes:
            return frame_class


def _set_wcs_aux_obs_coord(wcs, obs_frame):
    """
    Set (in-place) observer coordinate information on a WCS.

    Parameters
    ----------
    wcs : astropy.wcs.WCS
    obs_frame : astropy.coordinates.SkyCoord, astropy.coordinates.CoordinateFrame
    """
    # Sometimes obs_coord can be a SkyCoord, so convert down to a frame
    if hasattr(obs_frame, 'frame'):
        obs_frame = obs_frame.frame

    if isinstance(obs_frame, HeliographicStonyhurst):
        wcs.wcs.aux.hgln_obs = obs_frame.lon.to_value(u.deg)
    elif isinstance(obs_frame, HeliographicCarrington):
        wcs.wcs.aux.crln_obs = obs_frame.lon.to_value(u.deg)
    else:
        raise ValueError('obs_coord must be in a Stonyhurst or Carrington frame')
    # These two keywords are the same for Carrington and Stonyhurst
    wcs.wcs.aux.hglt_obs = obs_frame.lat.to_value(u.deg)
    wcs.wcs.aux.dsun_obs = obs_frame.radius.to_value(u.m)


def solar_frame_to_wcs_mapping(frame, projection='TAN'):
    """
    For a given frame, this function returns the corresponding WCS object.
    It registers the WCS coordinates types from their associated frame in the
    `astropy.wcs.utils.celestial_frame_to_wcs` registry.

    Parameters
    ----------
    frame : astropy.coordinates.BaseCoordinateFrame
    projection : str, optional

    Returns
    -------
    astropy.wcs.WCS
    """
    wcs = WCS(naxis=2)

    if hasattr(frame, 'rsun'):
        wcs.wcs.aux.rsun_ref = frame.rsun.to_value(u.m)

    if hasattr(frame, 'observer') and frame.observer is not None:
        if isinstance(frame.observer, BaseCoordinateFrame):
            observer = frame.observer
        elif frame.observer == 'self':
            observer = frame
        _set_wcs_aux_obs_coord(wcs, observer)

    if isinstance(frame, SunPyBaseCoordinateFrame):

        if frame.obstime:
            wcs.wcs.dateobs = frame.obstime.utc.isot

        if isinstance(frame, Helioprojective):
            xcoord = 'HPLN' + '-' + projection
            ycoord = 'HPLT' + '-' + projection
            wcs.wcs.cunit = ['arcsec', 'arcsec']
        elif isinstance(frame, Heliocentric):
            xcoord = 'SOLX'
            ycoord = 'SOLY'
            wcs.wcs.cunit = ['deg', 'deg']
        elif isinstance(frame, HeliographicCarrington):
            xcoord = 'CRLN' + '-' + projection
            ycoord = 'CRLT' + '-' + projection
            wcs.wcs.cunit = ['deg', 'deg']
        elif isinstance(frame, HeliographicStonyhurst):
            xcoord = 'HGLN' + '-' + projection
            ycoord = 'HGLT' + '-' + projection
            wcs.wcs.cunit = ['deg', 'deg']

    else:
        return None

    wcs.wcs.ctype = [xcoord, ycoord]

    return wcs


astropy.wcs.utils.WCS_FRAME_MAPPINGS.append([solar_wcs_frame_mapping])
astropy.wcs.utils.FRAME_WCS_MAPPINGS.append([solar_frame_to_wcs_mapping])
