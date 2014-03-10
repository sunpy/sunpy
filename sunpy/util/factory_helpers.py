# -*- coding: utf-8 -*-
"""
Some handy utils for factory classes
"""
import os
import datetime
import urllib2

from astropy.utils import isiterable

import sunpy.time

__all__ = ['is_url', 'is_time_string', 'is_time_date',
           'is_file', 'is_iterable_not_str']

def is_url(arg):
    try:
        urllib2.urlopen(arg)
    except ValueError:
        return False
    return True

def is_time_string(arg):
    try:
        sunpy.time.parse_time(arg)
    except (ValueError, TypeError):
        return False
    return True

def is_time_date(arg):
    """
    Checks if time is timerange time string or time strinf pair or datetime.date
    """
    if isiterable(arg) and len(arg) <= 2 and all([is_time_string(x) for x in arg]) or \
       isinstance(arg, sunpy.time.TimeRange) or isinstance(arg, datetime.date) \
       or is_time_string(arg):
           return True

    return False

def is_file(arg):
    return isinstance(arg,basestring) and os.path.isfile(os.path.expanduser(arg))

def is_iterable_not_str(arg):
    return isiterable(arg) and not isinstance(arg, basestring)