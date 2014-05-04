"""
CGS values of solar physics constants.
"""
from __future__ import absolute_import

from astropy.constants import Constant
import numpy as np
from . import _si

__all__ = ['physical_constants']

# This is really just to make the name shorter, so we can stick to a 
# maximum of 79 characters per line.
si_consts = _si.physical_constants

physical_constants = {}

# Not available in _si.py
#physical_constants['center density'] = (si_consts['center density'] * 
#                                        1.0e-3, 'g cm^-3', -1)
#physical_constants['center temperature'] = si_consts['center temperature']
#physical_constants['absolute magnitude'] = si_consts['absolute magnitude']
#physical_constants['mean energy production'] = (
#    si_consts['mean energy production'] * 1.0e4, 'erg g^-1', -1)
#physical_constants['ellipticity'] = si_consts['ellipticity']
#physical_constants['GM'] = (si_consts['GM'] * 1.0e6, 'cm^3 s^-2', -1)

# Not availabe in _si.py
#standard_model = _si.standard_model
