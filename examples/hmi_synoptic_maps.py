"""
HMI Daily Synoptic Map
----------------------

In this example we load the Daily Synoptic Maps produced by the HMI team. This
data is an interesting demonstration of SunPy's Map class as it is not in the
more common Helioprojective coordinate system, it is in Heliographic Carrington
coordinates and in a non-trivial Cylindrical Equal Area projection.

This example plots the HMI Daily Synoptic Maps, the file used in this example
can be downloaded from
`here <http://jsoc.stanford.edu/data/hmi/synoptic/hmi.Mldailysynframe_720s_nrt.fits>`_
"""

import numpy as np
import matplotlib.pyplot as plt

import astropy.units as u
from astropy.utils.data import download_file

import sunpy.map

###############################################################################
# Use astropy to download the file to a temp location.

filename = download_file('http://jsoc.stanford.edu/data/hmi/synoptic/hmi.Synoptic_Mr.2191.fits', cache=True)


###############################################################################
# We read this file in as a Map.

syn_map = sunpy.map.Map(filename)

###############################################################################
# There are a couple of oddities with this file, firstly the value of 'CUNIT2':

print(syn_map.meta['CUNIT2'])


###############################################################################
# That is not a unit! What this is telling us is that the latitude coordinate
# is actually the sine of latitude. According to the Thompson (2006) paper,
# CUNIT2 should be in degrees and CDELT2 should be multiplied by 180/pi. Also
# the value of CDELT1 has the wrong sign.

syn_map.meta['CUNIT2'] = 'degree'
syn_map.meta['CDELT2'] = 180/np.pi * syn_map.meta['CDELT2']
syn_map.meta['CDELT1'] *= -1

###############################################################################
# Now we create a SunPy Map from the data and header:

# Set the colorbar properties.
syn_map.plot_settings['cmap'] = 'hmimag'
syn_map.plot_settings['norm'] = plt.Normalize(-1500, 1500)


###############################################################################
# Create a figure with the Map's projection:
fig = plt.figure(figsize=(12, 5))
axes = plt.subplot(projection=syn_map)

# Plot the image
im = syn_map.plot()

# Set up the Sine Latitude Grid
x = axes.coords[0]
y = axes.coords[1]

x.set_coord_type('longitude', coord_wrap=360.)
y.set_coord_type('latitude')

x.set_major_formatter('dd')
y.set_major_formatter('dd')

x.set_axislabel("Carrington Longitude [deg]")
y.set_axislabel("Latitude [deg]")


x.set_ticks(color='black', exclude_overlapping=True)
y.set_ticks(color='black', exclude_overlapping=True)

# Hide the grid
axes.coords.grid(color='black', alpha=0.6, linestyle='dotted',
                 linewidth=0.5)

# Create a colorbar
cb = plt.colorbar(im, fraction=0.019, pad=0.1)
cb.set_label("LOS Magnetic Field [Gauss]")


# Another horrible hack to make the ticks draw on the RHS
axes.set_ylim((1, syn_map.data.shape[0]-1))

plt.title("{} {}-{}".format(syn_map.meta['content'], syn_map.meta['CAR_ROT'], syn_map.meta['CAR_ROT']+1))

plt.show()
