"""
===================================
Getting the location of the planets
===================================

How to get the position of planetary bodies im the solar system using
`astropy's solar system ephemeris <http://docs.astropy.org/en/stable/coordinates/solarsystem.html#solar-system-ephemerides>`__ information and sunpy.
"""
import matplotlib.pyplot as plt

from astropy.time import Time

from sunpy.coordinates import get_body_heliographic_stonyhurst

##############################################################################
# Lets grab the positions of each of the planets in Heliographic Stonyhurst
# coordinates.

obstime = Time('2014-05-15T07:54:00.005')
planet_list = ['earth', 'venus', 'mars', 'mercury', 'jupiter', 'neptune', 'uranus', 'sun']
planet_coord = [get_body_heliographic_stonyhurst(
    this_planet, time=obstime) for this_planet in planet_list]

##############################################################################
# Let's plot the results. Remember the Sun is at the center of this coordinate
# system.

fig = plt.figure()
ax = fig.add_subplot(projection='polar')
for this_planet, this_coord in zip(planet_list, planet_coord):
    ax.plot(this_coord.lon.to('rad'), this_coord.radius, 'o', label=this_planet)
ax.legend()

plt.show()
