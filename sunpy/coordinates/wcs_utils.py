# -*- coding: utf-8 -*-
from .frames import *

def solar_wcs_frame_mapping(wcs):
    xcoord = wcs.wcs.ctype[0][0:4]
    ycoord = wcs.wcs.ctype[1][0:4]

    if xcoord == 'HPLN' and ycoord == 'HPLT':
        return HelioProjective

    if xcoord == 'HGLN' and ycoord == 'HGLT':
        return HelioGraphicStonyhurst

    if xcoord == 'CRLN' and ycoord == 'CRLT':
        return HelioGraphicCarrington

    if xcoord == 'SOLX' and ycoord == 'SOLY':
        return HelioCentric

try:
    import astropy.wcs.utils

    astropy.wcs.utils.WCS_FRAME_MAPPINGS.append([solar_wcs_frame_mapping])

except ImportError:
    pass