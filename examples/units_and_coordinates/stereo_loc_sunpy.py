"""
===================================
Getting the location of STEREO-A using sunpy
===================================

How to get the position of planetary bodies im the solar system using
`astropy's solar system ephemeris <http://docs.astropy.org/en/stable/coordinates/solarsystem.html#solar-system-ephemerides>`__ information and sunpy.
"""
import matplotlib.pyplot as plt
from astropy.time import Time
import datetime
from sunpy.coordinates import get_body_heliographic_stonyhurst

today = datetime.datetime.now()
obstime = Time(today)

from sunpy.coordinates.ephemeris import get_horizons_coord
aia = get_horizons_coord('STEREO-A', obstime)

planet_list = ['mars', 'sun', 'earth']
planet_coord = [get_body_heliographic_stonyhurst(
    this_planet, time=obstime) for this_planet in planet_list]
planet_list.append('STEREO-A')
planet_coord.append(aia)

fig = plt.figure()
ax = fig.add_subplot(projection='polar')
for this_planet, this_coord in zip(planet_list, planet_coord):
    ax.plot(this_coord.lon.to('rad'), this_coord.radius, 'o', label=this_planet)
ax.legend()
plt.show()