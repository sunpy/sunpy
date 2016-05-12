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

# references
gsfc_fact = "http://nssdc.gsfc.nasa.gov/planetary/factsheet/sunfact.html"
allen = "Allen's Astrophysical Quantities 4th Ed."
asplund = "Asplund et al. 2006"
fivian = "Fivian et al. 2008"

physical_constants['mass'] = astrocon.M_sun
physical_constants['radius'] = astrocon.R_sun
physical_constants['luminosity'] = astrocon.L_sun
physical_constants['mean distance'] = astrocon.au

# following needs error estimate if appropriate
physical_constants['perihelion distance'] = Constant('perihelion',
                                                     "Perihelion Distance",
                                                     1.471e11, 'm', 0,
                                                     allen, system='si')

# following needs error estimate if appropriate
physical_constants['aphelion distance'] = Constant('aphelion',
                                                   "Aphelion Distance",
                                                   1.521e11, 'm', 0,
                                                   allen, system='si')

physical_constants['age'] = Constant('age', "Age of the Sun",
                                     4.6e9, 'year', 0.1e9,
                                     allen, system='si')

# A solar flux  (sfu) is traditional measure of solar radio flux.
physical_constants['solar flux unit'] = Constant('sfu', "Solar flux unit",
                                                 1e-22, 'W m**-2 Hz**-1', 0,
                                                 allen, system='si')

# following needs error estimate if appropriate
physical_constants['visual magnitude'] = Constant('V',
                                                  "Apparent visual magnitude",
                                                  -26.75, '', 0,
                                                  allen, system='si')

# The Sun as viewed from Earth
physical_constants['average angular size'] = Constant('theta', "Semidiameter",
                                                      959.63, 'arcsec', 0,
                                                      allen, system='si')

# following needs error estimate if appropriate
physical_constants['surface area'] = Constant('A', "Surface area",
                                              6.087e18, 'm**2', 0,
                                              allen, system='si')

# following needs error estimate if appropriate
physical_constants['average density'] = Constant('rho', "Mean density",
                                                 1409, 'kg m**-3', 0,
                                                 allen, system='si')

# following needs error estimate if appropriate
physical_constants['surface gravity'] = Constant('g', "Surface gravity",
                                                 274, 'm s**-1', 0,
                                                 allen, system='si')

# following needs error estimate if appropriate
physical_constants['moment of inertia'] = Constant('I', "Moment of inertia",
                                                   5.7e54, 'kg m**-2', 0,
                                                   allen, system='si')

# following needs error estimate if appropriate
physical_constants['volume'] = Constant('V', "Volume",
                                        1.4122e27, 'm**3', 0,
                                        allen, system='si')

# following needs error estimate if appropriate
physical_constants['escape velocity'] = Constant('v',
                                                 "Escape velocity at surface",
                                                 6.177e5, 'm s**-1', 0,
                                                 allen, system='si')

physical_constants['oblateness'] = Constant('', "oblateness",
                                            8.01, 'marcsec', 0.14,
                                            fivian, system='si')

# the following constants need references and error estimates if appropriate
physical_constants['metallicity'] = Constant('Z', "Metallicity",
                                             0.0122, '', 0.0,
                                             asplund, system='si')

sunspot_cycle_exp = "Average duration of sunspot cycle"
physical_constants['sunspot cycle'] = Constant('', sunspot_cycle_exp,
                                               11.4, 'year', 0,
                                               "", system='si')

physical_constants['average intensity'] = Constant('I', "Mean Intensity",
                                                   2.009e7, 'W m**-2 sr**-1',
                                                   0,
                                                   "", system='si')

effect_temp_exp = "Effective black-body temperature"
physical_constants['effective temperature'] = Constant('T', effect_temp_exp,
                                                       5778.0, 'K', 0,
                                                       "", system='si')

mass_conv_exp = "Mass conversion rate"
physical_constants['mass conversion rate'] = Constant('dm/dt', mass_conv_exp,
                                                      4300e6, 'kg s**-1', 0,
                                                      "", system='si')

# following needs error estimate if appropriate
physical_constants['center density'] = Constant('rho_center', "Center density",
                                                1.622e5, 'kg m**-3', 0,
                                                gsfc_fact, system='si')

# following needs error estimate if appropriate
cent_temp_exp = "Center temperature"
physical_constants['center temperature'] = Constant('T_center', cent_temp_exp,
                                                    1.571e7, 'K', 0,
                                                    gsfc_fact, system='si')

# following needs error estimate if appropriate
abs_magn_exp = "Absolute magnitude"
physical_constants['absolute magnitude'] = Constant('M_abs', abs_magn_exp,
                                                    +4.83, '', 0,
                                                    gsfc_fact, system='si')
# following needs error estimate if appropriate
mean_energy_exp = "mean energy production"
physical_constants['mean energy production'] = Constant('', mean_energy_exp,
                                                        193.7e-6, 'J kg**-1',
                                                        0,
                                                        gsfc_fact, system='si')

# following needs error estimate if appropriate
physical_constants['ellipticity'] = Constant('', "ellipticity",
                                             5e-5, '', 0,
                                             gsfc_fact, system='si')

# following needs error estimate if appropriate
physical_constants['GM'] = Constant('mu', "standard gravitational parameter",
                                    132.712e6, 'km**3 s**-2', 0,
                                    gsfc_fact, system='si')
