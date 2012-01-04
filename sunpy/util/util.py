"""Provides utility programs.

    Notes: 
    The astronomy-type utilities should probably be separated out into
    another file. 
    --schriste

"""

from __future__ import absolute_import
from scipy.constants import constants as con

__all__ = ["toggle_pylab", "anytim", "julian_day", "julian_centuries", 
           "day_of_year", "break_time", "degrees_to_hours", "degrees_to_arc",
           "kelvin_to_keV", "keV_to_kelvin", "unique"]

from matplotlib import pyplot
from datetime import datetime
from datetime import timedelta
import numpy as np
from itertools import izip, imap

# The number of days between Jan 1 1900 and the Julian reference date of 
# 12:00 noon Jan 1, 4713 BC
JULIAN_DAY_ON_NOON01JAN1900 = 2415021.0

def toggle_pylab(fn):
    """ A decorator to prevent functions from opening matplotlib windows
        unexpectedly when sunpy is run in interactive shells like ipython 
        --pylab. 

        Toggles the value of matplotlib.pyplot.isinteractive() to preserve the
        users' expections of pylab's behaviour in general. """

    if pyplot.isinteractive():
        def fn_itoggle(*args, **kwargs):
            pyplot.ioff()
            ret = fn(*args, **kwargs)
            pyplot.ion()
            return ret
        return fn_itoggle
    else:
        return fn

def anytim(time_string=None):
    """Given a time string will parse and return a datetime object.
    If no string is given then returns the datetime object for the current time.
    If a datetime object is passed in by mistake then it returns it without an 
    error.
    
    TODO: add ability to parse tai (International Atomic Time seconds since 
    Jan 1, 1958)
    """
    if time_string is None:
        return datetime.now()
    elif isinstance(time_string, datetime):
        return time_string
    elif isinstance(time_string, tuple):
        return datetime(*time_string)
    elif isinstance(time_string, int) or isinstance(time_string, float):
        return datetime(1979, 1, 1) + timedelta(0, time_string)
    else:
        time_format_list = \
            ["%Y-%m-%dT%H:%M:%S.%f",    # Example 2007-05-04T21:08:12.999999
             "%Y%m%dT%H%M%S.%f",        # Example 20070504T210812.999999
             "%Y/%m/%d %H:%M:%S",       # Example 2007/05/04 21:08:12
             "%Y/%m/%d %H:%M",          # Example 2007/05/04 21:08
             "%Y/%m/%d %H:%M:%S.%f",    # Example 2007/05/04 21:08:12.999999
             "%Y-%m-%d %H:%M:%S.%f",    # Example 2007-05-04 21:08:12.999999
             "%Y-%m-%dT%H:%M:%S.%fZ",   # Example 2007-05-04T21:08:12.999Z
             "%Y-%m-%d %H:%M:%S",       # Example 2007-05-04 21:08:12
             "%Y-%m-%dT%H:%M:%S",       # Example 2007-05-04T21:08:12
             "%Y-%m-%d %H:%M",          # Example 2007-05-04 21:08
             "%Y%m%dT%H%M%S",           # Example 20070504T210812
             "%Y-%b-%d %H:%M:%S",       # Example 2007-May-04 21:08:12
             "%Y-%b-%d %H:%M",          # Example 2007-May-04 21:08
             "%Y-%b-%d",                # Example 2007-May-04
             "%Y-%m-%d",                # Example 2007-05-04
             "%Y/%m/%d",                # Example 2007/05/04
             "%Y%m%d_%H%M%S"]           # Example 20070504_210812
        for time_format in time_format_list: 
            try: 
                return datetime.strptime(time_string, time_format)
            except ValueError:
                pass
    
        raise ValueError("%s is not a valid time string!" % time_string)

def julian_day(t=None):
    """Returns the (fractional) Julian day defined as the number of days 
    between the queried day and the reference date of 12:00 (noon) Jan 1, 4713 
    BC."""
    # Good online reference for fractional julian day
    # http://www.stevegs.com/jd_calc/jd_calc.htm
    
    JULIAN_REF_DAY = anytim('1900/1/1 12:00:00')
    time = anytim(t)
    
    tdiff = time - JULIAN_REF_DAY
    
    julian = tdiff.days + JULIAN_DAY_ON_NOON01JAN1900
   
    result = julian + 1/24.*(time.hour + time.minute/60.0 + 
                             time.second/(60.*60.))

    # This is because the days in datetime objects start at 00:00, 
    # not 12:00 as for Julian days.
    if time.hour >= 12:
        result = result - 0.5
    else:
        result = result + 0.5

    return result

def julian_centuries(t=None):
    """Returns the number of Julian centuries since 1900 January 0.5."""
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

def break_time(t=None):
    """Given a time returns a string. Useful for naming files."""
    #TODO: should be able to handle a time range
    time = anytim(t)
    result = time.strftime("%Y%m%d_%H%M%S")
    return result

def degrees_to_hours(angle):
    """Converts an angle from the degree notation to the hour, arcmin, arcsec 
    notation (returned as a tuple)."""
    hour = int(np.floor(angle / 15))
    remainder = angle / 15.0 - hour
    arcminute = int(np.floor(remainder * 60))
    remainder =  remainder * 60 - arcminute
    arcsecond = remainder * 60.0
    return [hour, arcminute, arcsecond]

def degrees_to_arc(angle):
    """Converts decimal degrees to degree, arcminute, 
    arcsecond (returned as a tuple)."""
    degree = int(np.floor(angle))
    remainder = angle - degree
    arcminute = int(np.floor(remainder * 60))
    remainder =  remainder * 60 - arcminute
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
    C = 299792458.
    ANGSTROM = units['Angstrom'][1]  
    try:
        type_, n = units[unit]
    except KeyError:
        raise ValueError('Cannot convert %s to Angstrom' % unit)
    
    if type_ == 'wavelength':
        x = n / ANGSTROM
        return value / x
    elif type_ == 'frequency':
        x = 1 / ANGSTROM / n
        return x * (C / value)
    elif type_ == 'energy':
        x = 1 / (ANGSTROM / 1e-2) / n
        return x * (1 / (8065.53 * value))
    else:
        raise ValueError('Unable to convert %s to Angstrom' % type_)

def kelvin_to_keV(temperature):
    """Convert from temperature expressed in Kelvin to a 
    temperature expressed in keV"""
    return temperature / (con.e / con.k * 1000.0) 

def keV_to_kelvin(temperature):
    """Convert from temperature expressed in keV to a temperature 
    expressed in Kelvin"""
    return temperature * (con.e / con.k * 1000.0) 

def unique(itr, key=None):
    items = set()
    if key is None:
        for elem in itr:
            if elem not in items:
                yield elem
                items.add(elem)
    else:
        for elem in itr:
            x = key(elem)
            if x not in items:
                yield elem
                items.add(x)

def print_table(lst, colsep=' ', linesep='\n'):
    width = [max(imap(len, col)) for col in izip(*lst)]
    return linesep.join(
        colsep.join(
            col.ljust(n) for n, col in izip(width, row)
        ) for row in lst
    )
