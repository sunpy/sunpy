# coding: utf-8
"""
===================================
Getting the location of the planets
===================================

How to get the position of planetary bodies im the solar system using
astropy ephemeris information and sunpy.
"""
import matplotlib.pyplot as plt
import numpy as np

from astropy.time import Time

from sunpy.coordinates import get_body_heliographic_stonyhurst

##############################################################################
# Lets grab the positions of each of the planets in stonyhurt coordinates.
obstime = Time('2014-05-15T07:54:00.005')
planet_list = ['earth', 'venus', 'mars', 'mercury', 'jupiter', 'neptune', 'uranus', 'sun']
planet_coord = [get_body_heliographic_stonyhurst(this_planet, time=obstime) for this_planet in planet_list]

##############################################################################
# Let's plot the results. Remember the Sun is at the center of this coordinate
# coordinate system.
ax = plt.subplot(1, 1, 1, projection='polar')
for this_planet, this_coord in zip(planet_list, planet_coord):
    plt.polar(np.deg2rad(this_coord.lon), this_coord.radius, 'o', label=this_planet)
plt.legend()
plt.show()
