# -*- coding: utf-8 -*-
import astropy.wcs.utils

from .frames import *

__all__ = ['solar_wcs_frame_mapping']

def solar_wcs_frame_mapping(wcs):
    """
    This function registers the coordinates frames to their FITS-WCS coordinate
    type values in the `astropy.wcs.utils.wcs_to_celestial_frame` registry.
    """

    xcoord = wcs.wcs.ctype[0][0:4]
    ycoord = wcs.wcs.ctype[1][0:4]

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

    if xcoord == 'HPLN' and ycoord == 'HPLT':
        return HelioProjective(dateobs=dateobs, L0=hglon, B0=hglat, D0=dsun)

    if xcoord == 'HGLN' and ycoord == 'HGLT':
        return HelioGraphicStonyhurst(dateobs=dateobs)

    if xcoord == 'CRLN' and ycoord == 'CRLT':
        return HelioGraphicCarrington(dateobs=dateobs)

    if xcoord == 'SOLX' and ycoord == 'SOLY':
        return HelioCentric(dateobs=dateobs, L0=hglon, B0=hglat, D0=dsun)

astropy.wcs.utils.WCS_FRAME_MAPPINGS.append([solar_wcs_frame_mapping])
