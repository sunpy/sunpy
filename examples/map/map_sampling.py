"""
=============================
Sampling a map at coordinates
=============================

How to sample a map at given coordinates.
"""
import matplotlib.pyplot as plt

import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.visualization import quantity_support

import sunpy.map
from sunpy.coordinates.utils import GreatArc
from sunpy.data.sample import AIA_171_IMAGE

quantity_support()
###############################################################################
# We start with the sample data.

m = sunpy.map.Map(AIA_171_IMAGE)

###############################################################################
# To set the coordinates, we'll use a great arc.

start = SkyCoord(735 * u.arcsec, -400 * u.arcsec, frame=m.coordinate_frame)
end = SkyCoord(735 * u.arcsec, -100 * u.arcsec, frame=m.coordinate_frame)
great_arc = GreatArc(start, end)
coords = great_arc.coordinates()

###############################################################################
# Plot the coordinates on the Sun.

fig = plt.figure()
ax = fig.add_subplot(projection=m)
m.plot(axes=ax, clip_interval=(1, 99.99)*u.percent)
ax.plot_coord(coords, color='c')


###############################################################################
# Plot the intensity along the arc from the start to the end point.
# This also shows how different interpolation methods can effect the result

fig, ax = plt.subplots(tight_layout=True)
for method in ['nearest', 'linear', 'cubic']:
    intensity = sunpy.map.sample_at_coords(m, coords, method=method)
    separation = coords.separation(coords[0]).to(u.arcsec)
    ax.plot(separation, intensity, label=method)

ax.set_xlabel(f'Separation from start of arc [{separation.unit}]')
ax.set_ylabel(f'Intensity [{intensity.unit}]')
ax.legend(title="Interpolation method")

plt.show()
