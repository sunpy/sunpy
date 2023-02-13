"""
===================================================
Obtaining a spacecraft trajectory from JPL Horizons
===================================================

This example shows how to obtain the trajectory of a spacecraft from JPL Horizons
and plot it relative to other bodies in the solar system.

JPL `Horizons <https://ssd.jpl.nasa.gov/horizons/>`__ can return the locations of
planets and minor bodies (e.g., asteroids) in the solar system, and it can also
return the location of a variety of major spacecraft.

You will need `astroquery <https://astroquery.readthedocs.io/>`__ installed.
"""

import matplotlib.pyplot as plt

import astropy.units as u

from sunpy.time import parse_time

from sunpy.coordinates import get_body_heliographic_stonyhurst, get_horizons_coord

##############################################################################
# We use :func:`~sunpy.coordinates.get_horizons_coord` to query JPL Horizons
# for the trajectory of Parker Solar Probe (PSP).  Let's request 50 days on
# either side of PSP's 14th closest approach to the Sun.

perihelion_14 = parse_time('2022-12-11 13:16')
psp = get_horizons_coord('Parker Solar Probe',
                         {'start': perihelion_14 - 50 * u.day,
                          'stop': perihelion_14 + 50 * u.day,
                          'step': '180m'})

##############################################################################
# Now to get the location of solar bodies, we can use
# `~sunpy.coordinates.get_body_heliographic_stonyhurst` which uses astropy to calculate
# the barycentric position:

planets = ['Earth', 'Sun']
planet_coords = [get_body_heliographic_stonyhurst(planet, time=now) for planet in planets]

##############################################################################
# For the purposes of plotting on a Matplotlib polar plot, we create a short
# convenience function to extract the necessary values in the appropriate units.

def coord_to_polar(coord):
    return coord.lon.to_value('rad'), coord.radius.to_value('AU')

##############################################################################
# Finally, we will create a polar plot.

fig = plt.figure()
ax = fig.add_subplot(projection='polar')
for coord in trajectory_coords:
    ax.plot(coord.lon.to('rad'), coord.radius, 'y.', markersize=2)
for planet, coord in zip(planets, planet_coords):
    ax.plot(coord.lon.to('rad'), coord.radius, 'o', label=planet)
ax.plot(coord.lon.to('rad'), coord.radius, label='PSP', color='#c8c825')
ax.legend(loc='upper center')

plt.show()

##############################################################################
# There are other tools that enable a similar style of figure.
# `solarmach <https://github.com/jgieseler/solarmach#usage>`__ is one such example.
