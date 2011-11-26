"""
CGS values of solar physics constants.
"""
from __future__ import absolute_import

import scipy.constants as _cd
import numpy as np
from . import _constants as _si

physical_constants = {}

physical_constants['mass'] = (_si.physical_constants['mass'][0] * 1000.0, 
                              'g', -1)
physical_constants['radius'] = (_si.physical_constants['radius'][0] * 100.0,
                                'cm', -1)
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
physical_constants['center density'] = (_si.physical_constants['center density']
                                        [0] * 1.0e-3, 'g cm^-3', -1)

