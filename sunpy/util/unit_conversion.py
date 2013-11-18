from __future__ import absolute_import

import numpy as np
from scipy.constants import constants
from astropy import units

__all__ = ['degrees_to_hours', 'degrees_to_arc', 'kelvin_to_keV', 
           'keV_to_kelvin', 'to_angstrom']
           
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

def to_angstrom(value, unit):
    """Given a value with a unit (given in a string), convert to angstroms"""
    value_quantity = value * units.Unit(unit)
    return value_quantity.to(units.angstrom, equivalencies=units.spectral()).value

def kelvin_to_keV(temperature):
    """Convert from temperature expressed in Kelvin to a 
    temperature expressed in keV"""
    return temperature / (constants.e / constants.k * 1000.0) 

def keV_to_kelvin(temperature):
    """Convert from temperature expressed in keV to a temperature 
    expressed in Kelvin"""
    return temperature * (constants.e / constants.k * 1000.0) 