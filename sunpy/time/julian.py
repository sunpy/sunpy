from __future__ import absolute_import

from sunpy.time import parse_time
from astropy.time import Time

__all__ = ['julian_day', 'julian_centuries']

# The number of days between Jan 1 1900 and the Julian reference date of
# 12:00 noon Jan 1, 4713 BC
JULIAN_DAY_ON_NOON01JAN1900 = 2415021.0
# Epoch J2000.0 = 2000 January 1.5 TD = JDE 2451 545.0
JULIAN_DAY_EPHEM_ON_NOON01JAN2000 = 2451545.0

def julian_day(t='now'):
    """Returns the (fractional) Julian day defined as the number of days
    between the queried day and the reference date of 12:00 (noon) Jan 1, 4713
    BC."""

    # Good online reference for fractional julian day
    # http://www.stevegs.com/jd_calc/jd_calc.htm

    return Time(parse_time(t), scale = 'utc').jd

def julian_centuries(t='now'):
    """Returns the number of Julian centuries since 2000 January 1.5 TD."""
    DAYS_IN_YEAR = 36525.0

    jde = Time(parse_time(t), scale='utc').tt.jd
    return (jde - JULIAN_DAY_EPHEM_ON_NOON01JAN2000)/ DAYS_IN_YEAR
