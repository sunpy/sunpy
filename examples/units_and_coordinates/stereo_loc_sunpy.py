"""
============================================
Getting the location of STEREO-A using sunpy
============================================

This example shows how to get and plot the position of planetary bodies within the solar system using
`astropy's solar system ephemeris <http://docs.astropy.org/en/stable/coordinates/solarsystem.html#solar-system-ephemerides>`__ information and sunpy.
"""
##############################################################################
# Import the required modules.

import datetime

import matplotlib.pyplot as plt

from astropy.time import Time

from sunpy.coordinates import get_body_heliographic_stonyhurst
from sunpy.coordinates.ephemeris import get_horizons_coord

##############################################################################
# You can make use of `~datetime.datetime.now` to pass current date and time.
# Now can use `get_horizons_coord()` to get the coordinates of STEREO-A.

today = datetime.datetime.now()
obstime = Time(today)
aia = get_horizons_coord('STEREO-A', obstime)

##############################################################################
# Now We are using `get_body_heliographic_stonyhurst()` to get the coordinates
# of other bodies.
# We have not used `get_body_heliographic_stonyhurst()` to find the coordinates
# of STEREO-A as it's position and velocity cannot be calculated with the
# 'builtin' ephemeris.

planets = ['mars', 'sun', 'earth']
planet_coords = [get_body_heliographic_stonyhurst(
    planet, time=obstime) for planet in planets]
planet_list.append('STEREO-A')
planet_coord.append(aia)

##############################################################################
# Finally, we will create a polar plot of these locations.

fig = plt.figure()
ax = fig.add_subplot(projection='polar')
for this_planet, this_coord in zip(planet_list, planet_coord):
    ax.plot(this_coord.lon.to('rad'), this_coord.radius, 'o', label=this_planet)
ax.legend()

plt.show()
