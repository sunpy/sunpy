# -*- coding: utf-8 -*-
#
# Author: Steven Christe <steven.d.christe@nasa.gov>
#
# <License info will go here...>
"""Provides utility programs.

    Notes: 
    The astronomy-type utilities should probably be separated out into
    another file. 
    --schriste

"""

from datetime import datetime
import numpy as np

def anytim(time_string=None):
    """Given a time string will parse and return a datetime object.
    If no string is given then returns the datetime object for the current time.
    If a datetime object is passed in by mistake then it returns it without an error.
    
    TODO: add ability to parse tai (International Atomic Time seconds since Jan 1, 1958)
    """
    if time_string is None:
        return datetime.now()
    if type(time_string) is type(datetime.now()):
        return time_string
    else:
        time_format_list = \
            ["%Y-%m-%dT%H:%M:%S.%f",    # Example 2007-05-04T21:08:12.1000000
             "%Y/%m/%d %H:%M:%S",       # Example 2007/05/04 21:08:12
             "%Y-%m-%d %H:%M:%S.%f",    # Example 2007/05/04 21:08:12.1000000
             "%Y-%m-%d %H:%M:%S",       # Example 2007-05-04 21:08:12
             "%Y-%b-%d %H:%M:%S"]       # Example 2007-May-04 21:08:12 
        for time_format in time_format_list: 
            try: 
                return datetime.strptime(time_string, time_format)
            except:
                pass
    
        raise ValueError("Not a valid time string!")

def julian_day(t=None):
    """Returns the (fractional) Julian day."""
    # The number of days between Jan 1 1900 and the Julian
    # reference date of 12:00 noon Jan 1, 4713 BC
    JULIAN_DAY_ON_NOON01JAN1900 = 2415020.5
    days = day_of_year(t)
    result = days + JULIAN_DAY_ON_NOON01JAN1900
    return result

def julian_centuries(t=None):
    """Returns the number of Julian centuries since 1900 January 0.5"""
    # The number of days between Jan 1 1900 and the Julian
    # reference date of 12:00 noon Jan 1, 4713 BC
    JULIAN_DAY_ON_NOON01JAN1900 = 2415020.5
    DAYS_IN_YEAR = 36525.0
    result = (julian_day(t) - JULIAN_DAY_ON_NOON01JAN1900) / DAYS_IN_YEAR
    return result

def day_of_year(t=None):
    """Returns the day of year."""
    SECONDS_IN_DAY = 60*60*24.0
    time = anytim(t)
    time_diff = anytim(t) - datetime(time.year, 1, 1, 0, 0, 0)
    result = time_diff.days + time_diff.seconds/SECONDS_IN_DAY
    return result

def degrees_to_hours(angle):
    """Converts an angle from the degree notation to the hour, arcmin, arcsec 
    notation (returned as a tuple)."""
    hour = int(np.floor(angle / 15))
    remainder = angle / 15.0 - hour
    arcminute = int(np.floor(remainder * 60))
    remainder =  remainder*60 - arcminute
    arcsecond = remainder * 60.0
    return [hour, arcminute, arcsecond]

def degrees_to_arc(angle):
    """Converts decimal degrees to degree, arcminute, 
    arcsecond (returned as a tuple)."""
    degree = int(np.floor(angle))
    remainder = angle - degree
    arcminute = int(np.floor(remainder * 60))
    remainder =  remainder*60 - arcminute
    arcsecond = remainder * 60.0
    return [degree, arcminute, arcsecond]