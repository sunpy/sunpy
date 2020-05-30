"""
=====================
Histograming map data
=====================

How to inspect the histogram of the data of a map.
"""
# sphinx_gallery_thumbnail_number = 2

import matplotlib.pyplot as plt
import numpy as np

import astropy.units as u
from astropy.coordinates import SkyCoord

import sunpy.map
from sunpy.data.sample import AIA_171_IMAGE

###############################################################################
# We start with the sample data and create a cutout.
aia = sunpy.map.Map(AIA_171_IMAGE)
bottom_left = SkyCoord(-300 * u.arcsec, 0 * u.arcsec, frame=aia.coordinate_frame)
top_right = SkyCoord(100 * u.arcsec, 400 * u.arcsec, frame=aia.coordinate_frame)
aia_smap = aia.submap(bottom_left, top_right=top_right)
aia_smap.plot()

###############################################################################
# The image of a `~sunpy.map.GenericMap` is always available in the data attribute.
# Map also provides shortcuts to the image minimum and maximum values.
# Let's create a histogram of the data in this submap.
num_bins = 50
bins = np.linspace(aia_smap.min(), aia_smap.max(), num_bins)
hist, bin_edges = np.histogram(aia_smap.data, bins=bins)

###############################################################################
# Let's plot the histogram as well as some standard values such as mean
# upper, and lower value and the one-sigma range.
plt.figure()
# Note that we have to use .ravel() here to avoid matplotlib interpreting each
# row in the array as a different dataset to histogram
plt.hist(aia_smap.data.ravel(), bins=bins, label='Histogram', histtype='step')
plt.xlabel('Intensity')
plt.axvline(aia_smap.min(), label='Data min={:.2f}'.format(aia_smap.min()), color='black')
plt.axvline(aia_smap.max(), label='Data max={:.2f}'.format(aia_smap.max()), color='black')
plt.axvline(aia_smap.data.mean(),
            label='mean={:.2f}'.format(aia_smap.data.mean()), color='green')
one_sigma = np.array([aia_smap.data.mean() - aia_smap.data.std(),
                      aia_smap.data.mean() + aia_smap.data.std()])
plt.axvspan(one_sigma[0], one_sigma[1], alpha=0.3, color='green',
            label='mean +/- std = [{:.2f}, {:.2f}]'.format(
            one_sigma[0], one_sigma[1]))
plt.axvline(one_sigma[0], color='green')
plt.axvline(one_sigma[1], color='red')
plt.yscale('log')
plt.legend(loc=9)

###############################################################################
# Finally let's overplot the one-sigma contours
fig = plt.figure()
fig.add_subplot(projection=aia_smap)
aia_smap.plot()
levels = one_sigma / aia_smap.max() * u.percent * 100
aia_smap.draw_contours(levels=levels, colors=['blue'])
plt.show()
