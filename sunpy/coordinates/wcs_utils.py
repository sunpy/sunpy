# -*- coding: utf-8 -*-
from __future__ import absolute_import, division
import astropy.wcs.utils

from astropy.wcs import WCSSUB_CELESTIAL

from .frames import *

__all__ = ['solar_wcs_frame_mapping']

def solar_wcs_frame_mapping(wcs):
    """
    This function registers the coordinates frames to their FITS-WCS coordinate
    type values in the `astropy.wcs.utils.wcs_to_celestial_frame` registry.
    """

    # SunPy Map adds some extra attributes to the WCS object.
    # We check for them here, and default to None.
    dateobs = wcs.wcs.dateobs if wcs.wcs.dateobs else None
    hglon = None
    hglat = None
    dsun = None

    if hasattr(wcs, 'heliographic_longitude'):
        hglon = wcs.heliographic_longitude

    if hasattr(wcs, 'heliographic_latitude'):
        hglat = wcs.heliographic_latitude

    if hasattr(wcs, 'dsun'):
        dsun = wcs.dsun


    # First we try the Celestial sub, which rectifies the order.
    # It will return any thing matching ??LN*, ??LT*
    wcss = wcs.sub([WCSSUB_CELESTIAL])

    # If the SUB works, use it.
    if wcss.naxis == 2:
        wcs = wcss

    xcoord = wcs.wcs.ctype[0][0:4]
    ycoord = wcs.wcs.ctype[1][0:4]

    if xcoord == 'HPLN' and ycoord == 'HPLT':
        return Helioprojective(dateobs=dateobs, L0=hglon, B0=hglat, D0=dsun)

    if xcoord == 'HGLN' and ycoord == 'HGLT':
        return HeliographicStonyhurst(dateobs=dateobs)

    if xcoord == 'CRLN' and ycoord == 'CRLT':
        return HeliographicCarrington(dateobs=dateobs)

    if xcoord == 'SOLX' and ycoord == 'SOLY':
        return Heliocentric(dateobs=dateobs, L0=hglon, B0=hglat, D0=dsun)


astropy.wcs.utils.WCS_FRAME_MAPPINGS.append([solar_wcs_frame_mapping])

# The following is a patch for wcsaxes 0.6 and lower:
try:
    import wcsaxes.wcs_utils
    if hasattr(wcsaxes.wcs_utils, 'WCS_FRAME_MAPPINGS'):
        wcsaxes.wcs_utils.WCS_FRAME_MAPPINGS.append([solar_wcs_frame_mapping])
except ImportError:
    pass

    # Now we try for heliocentric without the sub.
