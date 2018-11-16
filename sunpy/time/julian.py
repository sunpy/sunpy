from astropy.time import Time
from sunpy.time import parse_time

__all__ = ['julian_day', 'julian_centuries']


def julian_day(t='now'):
    """
    Wrap a UTC -> JD conversion from astropy.
    """
    return parse_time(t).jd


def julian_centuries(t='now'):
    """Returns the number of Julian centuries since J1900.0 (noon on 1900 January 0)."""
    DAYS_IN_JULIAN_CENTURY = 36525.0

    # J1900.0 is 2415021.0
    return (parse_time(t).jd - 2415020.0) / DAYS_IN_JULIAN_CENTURY
