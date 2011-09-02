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
             "%Y-%b-%d %H:%M:%S",       # Example 2007-May-04 21:08:12
             "%Y-%b-%d",                # Example 2007-May-04
             "%Y-%m-%d",                # Example 2007-05-03
             "%Y/%m/%d"]                # Example 2007/05/04
        for time_format in time_format_list: 
            try: 
                return datetime.strptime(time_string, time_format)
            except ValueError:
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

wavelength = [
    ('Angstrom', 1e-10),
    ('nm', 1e-9),
    ('micron', 1e-6),
    ('micrometer', 1e-6),
    ('mm', 1e-3),
    ('cm', 1e-2),
    ('m', 1e-6),
]
energy = [
    ('eV', 1),
    ('keV', 1e3),
    ('MeV', 1e6),
]
frequency = [
    ('Hz', 1),
    ('kHz', 1e3),
    ('MHz', 1e6),
    ('GHz', 1e9),
]
units = {}
for k, v in wavelength:
    units[k] = ('wavelength', v)
for k, v in energy:
    units[k] = ('energy', v)
for k, v in frequency:
    units[k] = ('frequency', v)

def to_angstrom(value, unit):
    C = 299792458
    ANGSTROM = units['Angstrom'][1]  
    try:
        type_, n = units[unit]
    except KeyError:
        raise ValueError('Cannot convert %s to Angstrom' % unit)
    
    if type_ == 'wavelength':
        k = n / ANGSTROM
        return value / k
    if type_ == 'frequency':
        k = 1 / ANGSTROM / n
        return k * (C / value)
    if type_ == 'energy':
        k = 1 / (ANGSTROM / 1e-2) / n
        return k * (1 / (8065.53 * value))
    else:
        raise ValueError('Unable to convert %s to Angstrom' % type_)

def unique(itr, key=None):
    items = set()
    if key is None:
        for elem in itr:
            if elem not in items:
                yield elem
                items.add(elem)
    else:
        for elem in itr:
            k = key(elem)
            if k not in items:
                yield elem
                items.add(k)
