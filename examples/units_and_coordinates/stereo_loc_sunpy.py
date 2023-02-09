"""
==============================================
Getting the trajectory of STEREO-A using sunpy
==============================================

This example shows how to get and plot the position of planetary bodies within the solar system using
`astropy's solar system ephemeris <http://docs.astropy.org/en/stable/coordinates/solarsystem.html#solar-system-ephemerides>`__ information and sunpy.
"""
import datetime

import matplotlib.pyplot as plt

from astropy.time import Time

from sunpy.coordinates import get_body_heliographic_stonyhurst
from sunpy.coordinates.ephemeris import get_horizons_coord

##############################################################################
# You can make use of `~datetime.datetime.now` to pass current date and time.
# Now can use `get_horizons_coord() <https://docs.sunpy.org/en/stable/generated/api/sunpy.coordinates.get_horizons_coord.html>`__
# to get the coordinates of STEREO-A.

today = datetime.datetime.now()
obstime = Time(today)
aia = get_horizons_coord('STEREO-A', obstime)
trajectory_coords = get_horizons_coord('STEREO-A', {'start': '2012-11-11', 'stop': '2013-05-11', 'step': '180m'})

##############################################################################
# Now We are using `get_body_heliographic_stonyhurst()` to get the coordinates
# of other bodies.
# We have not used `get_body_heliographic_stonyhurst()` to find the coordinates
# of STEREO-A as it's position and velocity cannot be calculated with the
# 'builtin' ephemeris.

planets = ['earth', 'sun']
planet_coords = [get_body_heliographic_stonyhurst(
    planet, time=obstime) for planet in planets]
planets.append('STEREO-A')
planet_coords.append(aia)

##############################################################################
# Now we will create a polar plot of these coordinates.

fig = plt.figure()
ax = fig.add_subplot(projection='polar')
for coord in trajectory_coords:
    ax.plot(coord.lon.to('deg'), coord.radius, 'y.', markersize=2)
for planet, coord in zip(planets, planet_coords):
    ax.plot(coord.lon.to('rad'), coord.radius, 'o', label=planet)

##############################################################################
# You have to seperately insert the label.

ax.plot(coord.lon.to('deg'), coord.radius, label='Trajectory of STEREO-A', color='yellow')
ax.legend(loc = 'upper right')

plt.show()
