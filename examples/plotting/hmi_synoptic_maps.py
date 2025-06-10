"""
===========================
Loading an HMI synoptic map
===========================

In this example we load a synoptic map produced by the HMI team. This
data is an interesting demonstration of sunpy's Map class as it is not in the
more common Helioprojective coordinate system, but in heliographic Carrington
coordinates and a cylindrical equal area (CEA) projection.
"""
import matplotlib.pyplot as plt

from astropy.utils.data import download_file

import sunpy.map

###############################################################################
# Download the file and read it into a Map.

filename = download_file(
    'http://jsoc.stanford.edu/data/hmi/synoptic/hmi.Synoptic_Mr.2191.fits', cache=True)
syn_map = sunpy.map.Map(filename)

###############################################################################
# Plot the results.

fig = plt.figure(figsize=(12, 5))
ax = plt.subplot(projection=syn_map)
im = syn_map.plot(axes=ax)

ax.coords[0].set_axislabel("Carrington Longitude [deg]")
ax.coords[1].set_axislabel("Latitude [deg]")

ax.coords.grid(color='black', alpha=0.6, linestyle='dotted', linewidth=0.5)

cb = plt.colorbar(im, fraction=0.019, pad=0.1)
cb.set_label(f"Radial magnetic field [{syn_map.unit}]")

# In order to make the x-axis ticks show, the bottom y-limit has to be adjusted slightly
ax.set_ylim(bottom=0)
ax.set_title(f"{syn_map.meta['content']},\n"
             f"Carrington rotation {syn_map.meta['CAR_ROT']}")

plt.show()
