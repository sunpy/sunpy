"""
============================================
Getting the trajectory of Parker Solar Probe
============================================

This example shows how to plot the trajectory of Parker Solar Probe within the solar system.
"""

import matplotlib.pyplot as plt

from astropy.time import Time

from sunpy.coordinates import get_body_heliographic_stonyhurst
from sunpy.coordinates.ephemeris import get_horizons_coord

##############################################################################
# You can make use of `~datetime.datetime.now` to pass current date and time.
# Now can use `get_horizons_coord() <https://docs.sunpy.org/en/stable/generated/api/sunpy.coordinates.get_horizons_coord.html>`__
# to get the coordinates of Parker Solar Probe.

now = Time(datetime.datetime.now())
trajectory_coords = get_horizons_coord('Parker Solar Probe', {'start': '2021-10-11', 'stop': '2022-01-12', 'step': '180m'})

##############################################################################
# Now We are using `get_body_heliographic_stonyhurst()` to get the coordinates
# of other bodies such as sun ans earth.
# We have not used `get_body_heliographic_stonyhurst()` to find the coordinates
# of Parker Solar Probe as it's position and velocity cannot be calculated with the
# 'builtin' ephemeris.

planets = ['EARTH', 'SUN']
planet_coords = [get_body_heliographic_stonyhurst(
    planet, time=now) for planet in planets]

##############################################################################
# Now we will create a polar plot of these coordinates.

fig = plt.figure()
ax = fig.add_subplot(projection='polar')
for coord in trajectory_coords:
    ax.plot(coord.lon.to('rad'), coord.radius, 'y.', markersize=2)
for planet, coord in zip(planets, planet_coords):
    ax.plot(coord.lon.to('rad'), coord.radius, 'o', label=planet)
ax.plot(coord.lon.to('rad'), coord.radius, label='Trajectory of Parker Solar Probe', color='#c8c825')
ax.legend(loc='upper right')

plt.show()
