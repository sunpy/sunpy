# -*- coding: utf-8 -*-
"""
=====================
Histograming map data
=====================

How to inspect the histogram of the data of a map.
"""
# sphinx_gallery_thumbnail_number = 2

import numpy as np
import matplotlib.pyplot as plt

import astropy.units as u

import sunpy.map
from astropy.coordinates import SkyCoord
from sunpy.data.sample import AIA_171_IMAGE

###############################################################################
# We start with the sample data and create a cutout.
aia = sunpy.map.Map(AIA_171_IMAGE)
bottom_left = SkyCoord(-300 * u.arcsec, 0 * u.arcsec, frame=aia.coordinate_frame)
top_right = SkyCoord(100 * u.arcsec, 400 * u.arcsec, frame=aia.coordinate_frame)
aia_smap = aia.submap(bottom_left, top_right)

###############################################################################
# The image of a `~sunpy.map.GenericMap` is always available in the data attribute.
# Map also provides shortcuts to the image minimum and maximum values.
# Let's create a histogram of the data in this submap.
num_bins = 50
bins = np.linspace(aia_smap.min(), aia_smap.max(), num_bins)
hist, bins = np.histogram(aia_smap.data, bins=bins)
width = 0.7 * (bins[1] - bins[0])
x = (bins[:-1] + bins[1:]) / 2

###############################################################################
# Let's plot the histogram as well as some standard values such as mean
# upper, and lower value and the one-sigma range.
plt.figure()
plt.bar(x, hist, align='center', width=width, label='Histogram')
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
plt.legend()
plt.show()

###############################################################################
# Finally let's overplot the one-sigma contours
fig = plt.figure()
fig.add_subplot(projection=aia_smap)
aia_smap.plot()
levels = one_sigma / aia_smap.max() * u.percent * 100
aia_smap.draw_contours(levels=levels, colors=['blue'])
plt.show()
