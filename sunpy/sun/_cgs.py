"""
CGS values of solar physics constants.
"""
from __future__ import absolute_import

import numpy as np
from . import _si

__all__ = ['physical_constants']

# This is really just to make the name shorter, so we can stick to a 
# maximum of 79 characters per line.
si_consts = _si.physical_constants

physical_constants = {}

physical_constants['mass'] = (si_consts['mass'][0] * 1.0e3, 'g', -1)
physical_constants['radius'] = (si_consts['radius'][0] * 1.0e2, 'cm', -1)
physical_constants['diameter'] = (physical_constants['radius'][0] * 2.0, 
                                  'cm', -1)
physical_constants['volume'] = (4 / 3. * np.pi *
                                physical_constants['radius'][0] ** 3, 'cm^3', 
                                -1)
physical_constants['surface area'] = (4 * np.pi *
                                      physical_constants['radius'][0] ** 2,
                                      'cm^2', -1)
physical_constants['average density'] = (physical_constants['mass'][0] / 
                                         physical_constants['volume'][0],
                                         'g cm^-3', -1)
physical_constants['center density'] = (si_consts['center density'][0] * 
                                        1.0e-3, 'g cm^-3', -1)
physical_constants['mean intensity'] = (si_consts['mean intensity'][0] * 
                                        1.0e3, 'erg cm^-2 sr^-1', -1)
physical_constants['effective temperature'] = si_consts['effective temperature']
physical_constants['center temperature'] = si_consts['center temperature']
physical_constants['luminosity'] = (si_consts['luminosity'][0] * 1.0e7, 
                                    'erg s^-1', -1)
physical_constants['absolute magnitude'] = si_consts['absolute magnitude']
physical_constants['visual magnitude'] = si_consts['visual magnitude']
physical_constants['mass conversion rate'] = (si_consts['mass conversion rate']
                                              [0] * 1.0e3, 'g s^-1', -1)
physical_constants['mean energy production'] = (
    si_consts['mean energy production'][0] * 1.0e4, 'erg g^-1', -1)
physical_constants['ellipticity'] = si_consts['ellipticity']
physical_constants['GM'] = (si_consts['GM'][0] * 1.0e6, 'cm^3 s^-2', -1)
physical_constants['surface gravity'] = (si_consts['surface gravity'][0] *
                                         1.0e-3, 'g cm^-3', -1)
physical_constants['escape velocity'] = (si_consts['escape velocity'][0] *
                                         1.0e5, 'cm s^-1', -1)
physical_constants['sunspot cycle'] = si_consts['sunspot cycle']
physical_constants['metallicity'] = si_consts['metallicity']

physical_constants['solar flux unit'] = (si_consts['solar flux unit'][0] *
                                         1.0e3, 'erg s^-1 cm^-2 Hz^-1', -1)
physical_constants['average_angular_size'] = si_consts['average_angular_size']

standard_model = _si.standard_model
