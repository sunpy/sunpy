from __future__ import absolute_import

from sunpy.time import parse_time

__all__ = ['julian_day', 'julian_centuries']

# The number of days between Jan 1 1900 and the Julian reference date of
# 12:00 noon Jan 1, 4713 BC
JULIAN_DAY_ON_NOON01JAN1900 = 2415021.0

def julian_day(t='now'):
    """Returns the (fractional) Julian day defined as the number of days
    between the queried day and the reference date of 12:00 (noon) Jan 1, 4713
    BC."""

    # Good online reference for fractional julian day
    # http://www.stevegs.com/jd_calc/jd_calc.htm
    JULIAN_REF_DAY = parse_time('1900/1/1 12:00:00')
    time = parse_time(t)

    tdiff = time - JULIAN_REF_DAY

    julian = tdiff.days + JULIAN_DAY_ON_NOON01JAN1900

    result = julian + 1 / 24. * (time.hour + time.minute / 60.0 +
                                 time.second / (60. * 60.))

    # This is because the days in datetime objects start at 00:00,
    # not 12:00 as for Julian days.
    if time.hour >= 12:
        result = result - 0.5
    else:
        result = result + 0.5

    return result

def julian_centuries(t='now'):
    """Returns the number of Julian centuries since 1900 January 0.5."""
    DAYS_IN_YEAR = 36525.0

    return (julian_day(t) - JULIAN_DAY_ON_NOON01JAN1900) / DAYS_IN_YEAR
