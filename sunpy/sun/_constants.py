"""
This module provides a non-comprehensive collection of solar physical constants.
"""
# TODO: Need better sources for some constants as well as error values.
import astropy.constants.astropyconst20 as astrocon
import astropy.units as u
from astropy.constants import Constant
from astropy.time import Time

__all__ = ['physical_constants']

physical_constants = {}

# references
gsfc_fact = "https://nssdc.gsfc.nasa.gov/planetary/factsheet/sunfact.html"
allen = "Allen's Astrophysical Quantities 4th Ed."
archinal = 'Archinal et al. 2018'
asplund = "Asplund et al. 2006"
fivian = "Fivian et al. 2008"
meeus = "Meeus 1998 Astronomical Algorithms 2nd Ed."

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
                                                 274, 'm s**-2', 0,
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

# longitude of the prime meridian (without light travel time to Earth
# and aberration effects) is 84.176 degrees eastward at J2000
physical_constants['W_0'] = Constant('W_0',
                                     'longitude of the prime meridian (epoch J2000.0)',
                                     84.176, 'deg',
                                     0, archinal,
                                     system='si')

# the definitional (fixed) rotation rate of the Sun relative to the stars (i.e., sidereal rotation rate)
physical_constants['sidereal rotation rate'] = Constant('', 'sidereal rotation rate',
                                                        14.1844, 'deg day**-1', 0,
                                                        archinal, system='si')

# time in Julian Days (TT) of the start of the first Carrington rotation
first_carrington_rotation = Time(2398167.4, format='jd', scale='tt')
physical_constants['first Carrington rotation (JD TT)'] = Constant('',
                                                                   'first Carrington '
                                                                   'rotation (JD TT)',
                                                                   first_carrington_rotation.tt.jd,
                                                                   'day', 0.1,
                                                                   meeus, system='si')

# length of the mean Carrington rotation as seen from Earth
# the rotation rate of the Sun appears to be slower by the rate at which the Earth orbits the Sun
period = 1 / (physical_constants['sidereal rotation rate'] / (360*u.deg) - 1 / u.yr)
physical_constants['mean synodic period'] = Constant('',
                                                     'mean synodic period',
                                                     period.to_value('day'), 'day', 0,
                                                     archinal, system='si')

# Sun's north pole is oriented RA=286.13 deg, Dec=63.87 deg in ICRS and HCRS
physical_constants['alpha_0'] \
    = Constant('alpha_0',
               'right ascension (RA) of the north pole (epoch J2000.0)',
               286.13, 'deg', 0, archinal, system='si')
physical_constants['delta_0'] \
    = Constant('delta_0',
               'declination of the north pole (epoch J2000.0)',
               63.87, 'deg', 0, archinal, system='si')
