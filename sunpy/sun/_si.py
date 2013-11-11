"""
Collection of solar physical constants.

The list is not meant to be comprehensive, but just a convenient list for 
everyday use.

.. todo:: Need better sources for some constants as well as error values.

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
physical_constants['average angular size'] = Constant('theta', "Semidiameter", 959.63, 
                                                 'arcsec', 0, "Allen's Astrophysical Quantities 4th Ed.", system='si')

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

physical_constants['oblateness'] = Constant('v', "oblateness", 8.01, 
                                                 'marcsec', 0.14, "Fivian et al. 2008", system='si')  

#the following constants need references and error estimates if appropriate
physical_constants['metallicity'] = Constant('Z', "Metallicity", 0.0122, '', 0.0, 
                                                 'Asplund et al. 2006', system='si')

physical_constants['sunspot cycle'] = Constant('v', "Average duration of sunspot cycle", 11.4, 
                                                 'year', 0, "", system='si')                                      

physical_constants['average intensity'] = Constant('I', "Mean Intensity", 2.009e7, 
                                                 'W m**-2 sr**-1', 0, "", system='si')                                      
        
physical_constants['effective temperature'] = Constant('T', "The effective black-body temperature of the Sun in Kelvin. ", 5778.0, 
                                                 'K', 0, "", system='si')

physical_constants['mass conversion rate'] = Constant('dm/dt', "Mass conversion rate", 4300e6, 
                                                 'kg s**-1', 0, "", system='si')