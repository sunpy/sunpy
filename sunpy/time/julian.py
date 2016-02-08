from __future__ import absolute_import

from astropy.time import Time
from sunpy.time import parse_time

__all__ = ['julian_day', 'julian_centuries']


def julian_day(t='now'):
    """
    Wrap a UTC -> JD conversion from astropy.
    """
    return Time(parse_time(t)).jd

# The number of days between Jan 1 1900 and the Julian reference date of
# 12:00 noon Jan 1, 4713 BC
JULIAN_DAY_ON_NOON01JAN1900 = 2415021.0


def julian_centuries(t='now'):
    """Returns the number of Julian centuries since 1900 January 0.5."""
    DAYS_IN_YEAR = 36525.0

    return (julian_day(t) - JULIAN_DAY_ON_NOON01JAN1900) / DAYS_IN_YEAR
