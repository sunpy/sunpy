# -*- coding: utf-8 -*-
from __future__ import absolute_import, division

import astropy.wcs.utils
from astropy.wcs import WCSSUB_CELESTIAL

from .frames import Helioprojective, Heliocentric, HeliographicStonyhurst, HeliographicCarrington

__all__ = ['solar_wcs_frame_mapping']


def solar_wcs_frame_mapping(wcs):
    """
    This function registers the coordinates frames to their FITS-WCS coordinate
    type values in the `astropy.wcs.utils.wcs_to_celestial_frame` registry.
    """

    # SunPy Map adds some extra attributes to the WCS object.
    # We check for them here, and default to None.
    dateobs = wcs.wcs.dateobs if wcs.wcs.dateobs else None

    if hasattr(wcs, 'heliographic_observer'):
        observer = wcs.heliographic_observer
    else:
        observer = None

    # First we try the Celestial sub, which rectifies the order.
    # It will return any thing matching ??LN*, ??LT*
    wcss = wcs.sub([WCSSUB_CELESTIAL])

    # If the SUB works, use it.
    if wcss.naxis == 2:
        wcs = wcss

    xcoord = wcs.wcs.ctype[0][0:4]
    ycoord = wcs.wcs.ctype[1][0:4]

    if xcoord == 'HPLN' and ycoord == 'HPLT':
        return Helioprojective(dateobs=dateobs, observer=observer)

    if xcoord == 'HGLN' and ycoord == 'HGLT':
        return HeliographicStonyhurst(dateobs=dateobs)

    if xcoord == 'CRLN' and ycoord == 'CRLT':
        return HeliographicCarrington(dateobs=dateobs)

    if xcoord == 'SOLX' and ycoord == 'SOLY':
        return Heliocentric(dateobs=dateobs, observer=observer)


astropy.wcs.utils.WCS_FRAME_MAPPINGS.append([solar_wcs_frame_mapping])
