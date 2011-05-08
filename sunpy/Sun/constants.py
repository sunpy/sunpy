"""
Collection of solar physical constants and conversion factors.

All constants are in SI units.

The list is not meant to be comprehensive, but just a convenient list for everyday use.
"""

"""
Written: Steven Christe (7-May-2011)
Modified: 

physical constants: imported from Review of Particle Physics 2010 (page 102), NASA Sun Fact Sheet
(http://nssdc.gsfc.nasa.gov/planetary/factsheet/sunfact.html), and 
Wikipedia (http://en.wikipedia.org/wiki/Sun)
Use at own risk:
TODO: References should be to published or standard sources (i.e. NOT websites)
"""

import scipy.constants as _cd

au = astronomical_unit = _cd.au

mass = 1.9884e30
equatorial_radius = radius = 6.9551e8
equatorial_diameter = 2*radius
volume = 1.412e18
surface_area = 6.0877e12
average_density = density = 1.408e3
center_density = 1.622e5
equatorial_surface_gravity = surface_gravity = 274
mean_intensity = intensity = 2.009e7
effective_temperature = 5778
center_temperature = 1.57e7
luminosity = 3.8427e26
spectral_classification = 'G2V'
absolute_magnitude = 4.83
visual_brightness = -26.74
# kg/s
mass_conversion_rate = 4300e6
# J/kg
mean_energy_production = 0.1937e-3
escape_velocity = 617.6e3
ellipticity = 0.00005
# km^3/s^2
GM = 132712e6
