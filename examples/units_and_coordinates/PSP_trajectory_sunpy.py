"""
============================================
Getting the trajectory of Parker Solar Probe
============================================

This example shows how to plot the trajectory of Parker Solar Probe within the solar system.
"""

import matplotlib.pyplot as plt

from sunpy.time import parse_time

from sunpy.coordinates import get_body_heliographic_stonyhurst
from sunpy.coordinates.ephemeris import get_horizons_coord

##############################################################################
##############################################################################
# By using `~sunpy.coordinates.ephemeris.get_horizons_coord`, we can query the 
# Jet Propulsion Laboratory `Horizons System <https://ssd.jpl.nasa.gov/horizons/>`__ 
# to get the location of any body they support, but here we will use it to get
# the coordinates of Parker Solar Probe.

now = parse_time("now")
trajectory_coords = get_horizons_coord('Parker Solar Probe', {'start': '2021-10-11', 'stop': '2022-01-12', 'step': '180m'})

##############################################################################
# Now to get the location of solar bodies, we can use
`~sunpy.coordinates.get_body_heliographic_stonyhurst` which uses astropy to calculate
# the barycentric position:

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
