"""
Collection of solar physical constants. Most constants are in SI units.

The list is not meant to be comprehensive, but just a convenient list for 
everyday use. All derived values (e.g. escape velocity, average density) 
are calculated. Use at own risk.

References:
    Review of Particle Physics 2010 (page 102)
    NASA Sun Fact Sheet 
      (http://nssdc.gsfc.nasa.gov/planetary/factsheet/sunfact.html)
    Wikipedia (http://en.wikipedia.org/wiki/Sun)

TODO: References should be to published or standard sources, NOT websites.
TODO: Add solar atmosphere values to standard model.
"""
from __future__ import absolute_import

import scipy.constants as _cd
import numpy as np

__all__ = ['physical_constants']

physical_constants = {}

# physical_constants[name] = (val, units, uncert)
physical_constants['mass'] = (1.9884e30, 'kg', -1)
physical_constants['radius'] = (6.95508e8, 'm', -1)
physical_constants['diameter'] = (physical_constants['radius'][0] * 2.0, 
                                  'm', -1)
physical_constants['volume'] = (4 / 3. * np.pi * 
                                physical_constants['radius'][0] ** 3, 'm^3', -1)
physical_constants['surface area'] = (4 * np.pi * 
                                      physical_constants['radius'][0] ** 2, 
                                      'm^2', -1)
physical_constants['average density'] = (physical_constants['mass'][0] / 
                                         physical_constants['volume'][0], 
                                         'kg m^-3', -1)
physical_constants['center density'] = (1.622e5, 'kg m^-3', -1)
physical_constants['mean intensity'] = (2.009e7, 'W m^-2 sr^-1', -1)
physical_constants['effective temperature'] = (5778.0, 'K', -1)
physical_constants['center temperature'] = (1.57e7, 'K', -1)
physical_constants['luminosity'] = (3.8427e26, 'J s^-1', -1)
physical_constants['absolute magnitude'] = (4.83, 'None', -1)
physical_constants['visual magnitude'] = (-26.74, 'None', -1)
physical_constants['mass conversion rate'] = (4300e6, 'kg s^-1', -1)
physical_constants['mean energy production'] = (0.1937, 'J kg^-1', -1)
physical_constants['ellipticity'] = (0.00005, 'None', -1)
physical_constants['GM'] = (_cd.G * physical_constants['mass'][0], 'm^3 s^-2', 
                            -1)
physical_constants['surface gravity'] = (physical_constants['GM'][0] / 
                                         physical_constants['radius'][0] ** 2, 
                                         'kg m^-3', -1)
physical_constants['escape velocity'] = (1e-3 * np.sqrt(2 * 
                                         physical_constants['GM'][0] / 
                                         physical_constants['radius'][0]), 
                                         'km s^-1', -1)
physical_constants['sunspot cycle'] = (11.4, 'years', -1)
physical_constants['metallicity'] = (0.0122, 'None', -1)

# A solar flux unit (sfu) is traditional measure of solar radio flux.
physical_constants['solar flux unit'] = (1e-22, 'W m^-2 Hz^-1', 0)

# Solar radius measured outside earth's atmosphere in arcseconds
physical_constants['average_angular_size'] = (961.07064, 'arcsec', -1)

# Standard Model - Interior Structure
# adapted from Turck-Chieze et al. (1988)
# Composition X = 0.7046, Y = 0.2757, Z = 0.0197
standard_model = {}
standard_model['interior'] = {}
# Radius - unit R_sun
standard_model['interior']['radius'] = [0, 0.01, 0.022, 0.061, 0.090, 0.120,
                                        0.166, 0.202, 0.246, 0.281, 0.317,
                                        0.370, 0.453, 0.611, 0.7304, 0.862,
                                        0.965, 1.0000]
# mass - unit M_sun
standard_model['interior']['mass'] = [0, 0.0001, 0.002, 0.020, 0.057, 0.115,
                                      0.235, 0.341, 0.470, 0.562, 0.647, 0.748,
                                      0.854, 0.951, 0.9809, 0.9964, 0.9999,
                                      1.000]
# luminosity - unit L_sun
standard_model['interior']['luminosity'] = [0, 0.0009, 0.009, 0.154, 0.365, 
                                            0.594, 0.845, 0.940, 0.985, 0.997,
                                            0.992, 0.9996, 1.000, 1.0, 1.0, 1.0,
                                            1.0, 1.0]
# temperature - unit 10^6 K
standard_model['interior']['temperature'] = [15.513, 15.48, 15.36, 14.404, 
                                             13.37, 12.25, 10.53, 9.30, 8.035,
                                             7.214, 6.461, 5.531, 4.426, 2.981,
                                             2.035, 0.884, 0.1818, 0.005770]
# density - unit g cm^-3
standard_model['interior']['density'] = [147.74, 146.66, 142.73, 116.10, 93.35,
                                         72.73, 48.19, 34.28, 21.958, 15.157,
                                         10.157, 5.566, 2.259, 0.4483, 0.1528,
                                         0.042, 0.00361, 1.99e-7]

standard_model['evolution'] = {}
# time - unit 10^9 years
standard_model['evolution']['time'] = [0, 0.143, 0.856, 1.863, 2.193, 3.020,
                                       3.977, 4.587, 5.506, 6.074, 6.577, 7.027,
                                       7.728, 8.258, 8.7566, 9.805]
# luminosity - unit L_sun
standard_model['evolution']['luminosity'] = [0.7688, 0.7248, 0.7621, 0.8156,
                                             0.8352, 0.8855, 0.9522, 1.0, 1.079,
                                             1.133, 1.186, 1.238, 1.318, 1.399,
                                             1.494, 1.760]
# radius - unit R_sun
standard_model['evolution']['radius'] = [0.872, 0.885, 0.902, 0.924, 0.932,
                                         0.953, 0.981, 1.0, 1.035, 1.059, 1.082,
                                         1.105, 1.143, 1.180, 1.224, 1.361]
# central temperature - unit 10^6 K
standard_model['evolution']['temperature central'] = [13.35, 13.46, 13.68, 
                                                      14.08, 14.22, 14.60,
                                                      15.12, 15.51, 16.18,
                                                      16.65, 17.13, 17.62,
                                                      18.42, 18.74, 18.81,
                                                      19.25]
