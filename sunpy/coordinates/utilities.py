# -*- coding: utf-8 -*-
from __future__ import absolute_import, division
import datetime

from astropy.time import Time
from astropy.coordinates import SkyCoord, get_body_barycentric

from sunpy.time import parse_time

from .frames import HeliographicStonyhurst

__all__ = ['get_earth']


def get_earth(time='now'):
    """
    Return a SkyCoord for the location of the Earth at a specified time in the
    HeliographicStonyhurst frame.

    Parameters
    ----------
    time : various
        Time to use in a parse_time-compatible format

    Returns
    -------
    out : SkyCoord
        SkyCoord for the location of the Earth in the HeliographicStonyhurst frame
    
    """
    obstime = Time(parse_time(time))

    return SkyCoord(get_body_barycentric('earth', obstime), obstime=obstime).transform_to(HeliographicStonyhurst)
