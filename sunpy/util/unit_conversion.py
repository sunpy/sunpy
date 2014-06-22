from __future__ import absolute_import

import numpy as np
from scipy.constants import constants
from astropy import units

__all__ = ['kelvin_to_keV', 'keV_to_kelvin', 'to_angstrom']

boltz_unit = units.J / units.K

def to_angstrom(value, unit):
    """Given a value with a unit (given in a string), convert to angstroms"""
    value_quantity = value * units.Unit(unit)
    return value_quantity.to(units.angstrom, equivalencies=units.spectral()).value

def kelvin_to_keV(temperature):
    """Convert from temperature expressed in Kelvin to a 
    temperature expressed in keV."""
    if not isinstance(temperature, units.Quantity):
        raise ValueError("temperature should be kelvin Quantity")
    return temperature / ((constants.e * units.J / units.eV) / (constants.k * boltz_unit) * 1000.0) 
	
def keV_to_kelvin(temperature):
    """Convert from temperature expressed in keV to a temperature 
    expressed in Kelvin."""
    if not isinstance(temperature, units.Quantity):
        raise ValueError("temperature should be keV Quantity")
    return temperature * ((constants.e * units.J / units.eV) / (constants.k * boltz_unit) * 1000.0) 
