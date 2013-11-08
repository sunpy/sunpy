"""
Collection of solar physical constants. Most constants are in SI s.

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

from astropy.constants import Constant
import astropy.constants as astrocon

__all__ = ['physical_constants']

physical_constants = {}

physical_constants['mass'] = astrocon.M_sun
physical_constants['radius'] = astrocon.R_sun
physical_constants['luminosity'] = astrocon.L_sun
physical_constants['mean distance'] = astrocon.au

# following needs error estimate if appropriate
physical_constants['perihelion distance'] = Constant('perihelion', "Perihelion Distance", 1.471e11, 'm', 0, 
                                     "Allen's Astrophysical Quantities 4th Ed.", system='si')

# following needs error estimate if appropriate
physical_constants['aphelion distance'] = Constant('aphelion', "Aphelion Distance", 1.521e11, 'm', 0, 
                                     "Allen's Astrophysical Quantities 4th Ed.", system='si')

physical_constants['age'] = Constant('age', "Age of the Sun", 4.6e9, 'year', 0.1e9, 
                                     "Allen's Astrophysical Quantities 4th Ed.", system='si')

# A solar flux  (sfu) is traditional measure of solar radio flux.
physical_constants['solar flux unit'] = Constant('sfu', "Solar flux unit", 1e-22, 
                                                 'W m**-2 Hz**-1', 0, "Allen's Astrophysical Quantities 4th Ed.", system='si')

# following needs error estimate if appropriate
physical_constants['visual magnitude'] = Constant('V', "Apparent visual magnitude", -26.75, 
                                                 '', 0, "Allen's Astrophysical Quantities 4th Ed.", system='si')

# The Sun as viewed from Earth
physical_constants['average_angular_size'] = Constant('theta', "Semidiameter", 959.63, 
                                                 'arcsec', 0, "Allen's Astrophysical Quantities 4th Ed.", system='si')

# following needs error estimate if appropriate
physical_constants['surface area'] = Constant('A', "Surface area", 6.087e18, 
                                                 'm**2', 0, "Allen's Astrophysical Quantities 4th Ed.", system='si')

# following needs error estimate if appropriate
physical_constants['surface area'] = Constant('A', "Surface area", 6.087e18, 
                                                 'm**2', 0, "Allen's Astrophysical Quantities 4th Ed.", system='si')

# following needs error estimate if appropriate
physical_constants['average density'] = Constant('rho', "Mean density", 1409, 
                                                 'kg m**-3', 0, "Allen's Astrophysical Quantities 4th Ed.", system='si')

# following needs error estimate if appropriate
physical_constants['surface gravity'] = Constant('g', "Surface gravity", 274, 
                                                 'm s**-1', 0, "Allen's Astrophysical Quantities 4th Ed.", system='si')

# following needs error estimate if appropriate
physical_constants['moment of inertia'] = Constant('I', "Moment of inertia", 5.7e54, 
                                                 'kg m**-2', 0, "Allen's Astrophysical Quantities 4th Ed.", system='si')

# following needs error estimate if appropriate
physical_constants['volume'] = Constant('V', "Volume", 1.4122e27, 
                                                 'm**3', 0, "Allen's Astrophysical Quantities 4th Ed.", system='si')

# following needs error estimate if appropriate
physical_constants['escape velocity'] = Constant('v', "Escape velocity at surface", 6.177e5, 
                                                 'm s**-1', 0, "Allen's Astrophysical Quantities 4th Ed.", system='si')

#the following constants need references and error estimates if appropriate
physical_constants['metallicity'] = Constant('v', "Metallicity", 0.0122, 
                                                 '', 0, "", system='si')

physical_constants['sunspot cycle'] = Constant('v', "Average duration of sunspot cycle", 11.4, 
                                                 'year', 0, "", system='si')                                      

physical_constants['ellipticity'] = Constant('v', "Ellipticity", 0.00005, 
                                                 '', 0, "", system='si')  
physical_constants['average intensity'] = Constant('I', "Mean Intensity", 2.009e7, 
                                                 'W m**-2 sr**-1', 0, "", system='si')                                      
        
physical_constants['effective temperature'] = Constant('T', "Effective surface temperature", 5778.0, 
                                                 'K', 0, "", system='si')

physical_constants['mass conversion rate'] = Constant('dm/dt', "Mass conversion rate", 4300e6, 
                                                 'kg s**-1', 0, "", system='si')

# Solar radius measured outside earth's atmosphere in arcseconds

# Standard Model - Interior Structure
# adapted from Turck-Chieze et al. (1988)
# Composition X = 0.7046, Y = 0.2757, Z = 0.0197
standard_model = {}
standard_model['interior'] = {}
# Radius -  R_sun
standard_model['interior']['radius'] = [0, 0.01, 0.022, 0.061, 0.090, 0.120,
                                        0.166, 0.202, 0.246, 0.281, 0.317,
                                        0.370, 0.453, 0.611, 0.7304, 0.862,
                                        0.965, 1.0000]
# mass -  M_sun
standard_model['interior']['mass'] = [0, 0.0001, 0.002, 0.020, 0.057, 0.115,
                                      0.235, 0.341, 0.470, 0.562, 0.647, 0.748,
                                      0.854, 0.951, 0.9809, 0.9964, 0.9999,
                                      1.000]
# luminosity -  L_sun
standard_model['interior']['luminosity'] = [0, 0.0009, 0.009, 0.154, 0.365, 
                                            0.594, 0.845, 0.940, 0.985, 0.997,
                                            0.992, 0.9996, 1.000, 1.0, 1.0, 1.0,
                                            1.0, 1.0]
# temperature -  10^6 K
standard_model['interior']['temperature'] = [15.513, 15.48, 15.36, 14.404, 
                                             13.37, 12.25, 10.53, 9.30, 8.035,
                                             7.214, 6.461, 5.531, 4.426, 2.981,
                                             2.035, 0.884, 0.1818, 0.005770]
# density -  g cm^-3
standard_model['interior']['density'] = [147.74, 146.66, 142.73, 116.10, 93.35,
                                         72.73, 48.19, 34.28, 21.958, 15.157,
                                         10.157, 5.566, 2.259, 0.4483, 0.1528,
                                         0.042, 0.00361, 1.99e-7]

standard_model['evolution'] = {}
# time -  10^9 years
standard_model['evolution']['time'] = [0, 0.143, 0.856, 1.863, 2.193, 3.020,
                                       3.977, 4.587, 5.506, 6.074, 6.577, 7.027,
                                       7.728, 8.258, 8.7566, 9.805]
# luminosity -  L_sun
standard_model['evolution']['luminosity'] = [0.7688, 0.7248, 0.7621, 0.8156,
                                             0.8352, 0.8855, 0.9522, 1.0, 1.079,
                                             1.133, 1.186, 1.238, 1.318, 1.399,
                                             1.494, 1.760]
# radius -  R_sun
standard_model['evolution']['radius'] = [0.872, 0.885, 0.902, 0.924, 0.932,
                                         0.953, 0.981, 1.0, 1.035, 1.059, 1.082,
                                         1.105, 1.143, 1.180, 1.224, 1.361]
# central temperature -  10^6 K
standard_model['evolution']['temperature central'] = [13.35, 13.46, 13.68, 
                                                      14.08, 14.22, 14.60,
                                                      15.12, 15.51, 16.18,
                                                      16.65, 17.13, 17.62,
                                                      18.42, 18.74, 18.81,
                                                      19.25]
