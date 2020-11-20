import warnings

import astropy.units as u
import astropy.wcs.utils
from astropy.coordinates import BaseCoordinateFrame
from astropy.wcs import WCS

from sunpy.util.exceptions import SunpyUserWarning
from .frames import (
    BaseCoordinateFrame,
    Heliocentric,
    HeliographicCarrington,
    HeliographicStonyhurst,
    Helioprojective,
    SunPyBaseCoordinateFrame,
)

__all__ = ['solar_wcs_frame_mapping', 'solar_frame_to_wcs_mapping']


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

    observer = None
    for frame, attr_names in required_attrs.items():
        attrs = [getattr(wcs.wcs.aux, attr_name) for attr_name in attr_names]
        if all([attr is not None for attr in attrs]):
            kwargs = {'obstime': dateobs}
            if issubclass(frame, HeliographicCarrington):
                kwargs['observer'] = 'self'

            observer = frame(attrs[0] * u.deg,
                             attrs[1] * u.deg,
                             attrs[2] * u.m,
                             **kwargs)

    # Get rsun from the WCS auxillary information
    rsun = wcs.wcs.aux.rsun_ref
    if rsun is not None:
        rsun *= u.m

    # These custom attributes were always used in sunpy < 2.1; these warnings
    # can be converted into errors in sunpy 3.1
    if hasattr(wcs, 'heliographic_observer'):
        warnings.warn('Support for the .heliographic_observer attribute on a WCS is deprecated. '
                      'Set observer keywords in the FITS header, or directly set the wcs.wcs.aux '
                      'values instead.',
                      SunpyUserWarning)
        if observer is None:
            observer = wcs.heliographic_observer
        else:
            warnings.warn('Observer information present in WCS auxillary information, ignoring '
                          '.heliographic_observer')

    if hasattr(wcs, 'rsun'):
        warnings.warn('Support for the .rsun attribute on a WCS is deprecated. '
                      'Set observer keywords in the FITS header, or directly set the wcs.wcs.aux '
                      'values instead.',
                      SunpyUserWarning)
        if rsun is None:
            rsun = wcs.rsun
        else:
            warnings.warn('rsun information present in WCS auxillary information, ignoring '
                          '.rsun')

    # Truncate the ctype to the first four letters
    ctypes = {c[:4] for c in wcs.wcs.ctype}

    if {'HPLN', 'HPLT'} <= ctypes:
        return Helioprojective(obstime=dateobs, observer=observer, rsun=rsun)

    if {'HGLN', 'HGLT'} <= ctypes:
        return HeliographicStonyhurst(obstime=dateobs)

    if {'CRLN', 'CRLT'} <= ctypes:
        return HeliographicCarrington(obstime=dateobs, observer=observer)

    if {'SOLX', 'SOLY'} <= ctypes:
        return Heliocentric(obstime=dateobs, observer=observer)


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
    frame : astropy.coordiantes.BaseCoordinateFrame
    projection : str, optional

    Returns
    -------
    astropy.wcs.WCS
    """
    wcs = WCS(naxis=2)

    if hasattr(frame, 'rsun'):
        wcs.wcs.aux.rsun_ref = frame.rsun.to_value(u.m)

    if hasattr(frame, 'observer') and isinstance(frame.observer, BaseCoordinateFrame):
        _set_wcs_aux_obs_coord(wcs, frame.observer)

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
