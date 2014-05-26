from __future__ import absolute_import

import numpy as np
from scipy.constants import constants
from astropy import units

__all__ = ['degrees_to_hours', 'degrees_to_arc', 'kelvin_to_keV', 
           'keV_to_kelvin', 'to_angstrom']

boltz_unit = units.J / units.K
           
def degrees_to_hours(angle):
    """Converts an angle from the degree notation to the hour, arcmin, arcsec 
    notation (returned as a tuple)."""
    if not isinstance(angle, units.Quantity):
		raise ValueError("angle should be a Quantity")
    hour = (np.floor(angle / 15)) / units.deg
    remainder = ((angle / 15.0) / units.deg) - hour
    arcminute = (np.floor(remainder * 60)) 
    remainder =  remainder * 60 - arcminute
    arcsecond = remainder * 60.0
    return [hour * units.hourangle, arcminute * units.arcmin, arcsecond * units.arcsec]

def degrees_to_arc(angle):
    """Converts decimal degrees to degree, arcminute, 
    arcsecond (returned as a tuple)."""
    if not isinstance(angle, units.Quantity):
        raise ValueError("angle should be a Quantity")	
    degree = (np.floor(angle)) / units.deg
    remainder = (angle / units.deg) - degree
    arcminute = (np.floor(remainder * 60))
    remainder =  remainder * 60 - arcminute
    arcsecond = remainder * 60.0
    return [degree * units.degree, arcminute * units.arcmin, arcsecond * units.arcsec]

def to_angstrom(value, unit):
    """Given a value with a unit (given in a string), convert to angstroms"""
    value_quantity = value * units.Unit(unit)
    return value_quantity.to(units.angstrom, equivalencies=units.spectral()).value

def kelvin_to_keV(temperature):  #units are not converted properly
    """Convert from temperature expressed in Kelvin to a 
    temperature expressed in keV"""
	if not isinstance(temperature, units.Quantity):
	    raise ValueError("temperature should be kelvin Quantity")
    return temperature / ((constants.e * units.J) / (constants.k * boltz_unit) * 1000.0) 

def keV_to_kelvin(temperature):  #units are not converted properly
    """Convert from temperature expressed in keV to a temperature 
    expressed in Kelvin"""
	if not isinstance(temperature, units.Quantity):
	    raise ValueError("temperature should be eV Quantity")
    return temperature * ((constants.e * units.J) / (constants.k * units.J / units.K) * 1000.0) 
