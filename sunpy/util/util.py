"""Provides utility programs.

    Notes: 
    The astronomy-type utilities should probably be separated out into
    another file. 
    --schriste

"""

from __future__ import absolute_import
from scipy.constants import constants as con

__all__ = ["toggle_pylab", "degrees_to_hours", "degrees_to_arc",
           "kelvin_to_keV", "keV_to_kelvin", "unique", "print_table", 
           "to_angstrom"]

from matplotlib import pyplot
import numpy as np
from itertools import izip, imap

def to_signed(dtype):
    """ Return dtype that can hold data of passed dtype but is signed.
    Raise ValueError if no such dtype exists.
    
    Parameters
    ----------
    dtype : np.dtype
        dtype whose values the new dtype needs to be able to represent.
    """
    if dtype.kind == "u":
        if dtype.itemsize == 8:
            raise ValueError("Cannot losslessy convert uint64 to int.")
        dtype = "int%d" % (min(dtype.itemsize * 2 * 8, 64))
    return np.dtype(dtype)

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
