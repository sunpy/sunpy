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

    dateobs = wcs.wcs.dateobs or None

    # SunPy Map adds 'heliographic_observer' and 'rsun' attributes to the WCS
    # object. We check for them here, and default to None.
    if hasattr(wcs, 'heliographic_observer'):
        observer = wcs.heliographic_observer
    else:
        observer = None

    if hasattr(wcs, 'rsun'):
        rsun = wcs.rsun
    else:
        rsun = None

    # Truncate the ctype to the first four letters
    ctypes = {c[:4] for c in wcs.wcs.ctype}

    if {'HPLN', 'HPLT'} <= ctypes:
        return Helioprojective(obstime=dateobs, observer=observer, rsun=rsun)

    if {'HGLN', 'HGLT'} <= ctypes:
        return HeliographicStonyhurst(obstime=dateobs)

    if {'CRLN', 'CRLT'} <= ctypes:
        return HeliographicCarrington(obstime=dateobs)

    if {'SOLX', 'SOLY'} <= ctypes:
        return Heliocentric(obstime=dateobs, observer=observer)


astropy.wcs.utils.WCS_FRAME_MAPPINGS.append([solar_wcs_frame_mapping])
