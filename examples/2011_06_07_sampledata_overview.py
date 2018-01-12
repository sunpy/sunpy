"""
========================
Sample Data Set Overview
========================

An overview of the coordinated sample data set for the SOL-2011-06-07 event.
"""
from __future__ import print_function, division

import matplotlib.pyplot as plt

import sunpy.map
import sunpy.timeseries
import sunpy.data.sample as sample_data

###############################################################################
# On 2011 June 7, various solar instruments observed a spectacular solar
# eruptive event coming from NOAA AR 11226. The event included an M2.5 flare, a
# filament eruption, a CME, and a global coronal EUV wave (IAU standard:
# SOL2011-06-07T06:24:00L045C112). This event was spectacular because it
# features the ejection of a large amount of prominence material, much of which
# fails to escape and falls back to the solar surface.
# This example shows off our coordinated sample data set.
# This event got a little bit of press coverage (e.g. `National Geographics
# <http://news.nationalgeographic.com/news/2011/06/110608-solar-flare-sun-science-space/>`_,
# `Discover Magazine <http://blogs.discovermagazine.com/badastronomy/2011/06/07/the-sun-lets-loose-a-huge-explosion/#.WlVnNnWQVxB>`_)
# and the literature contains a number of a papers about it (e.g. `Li et al.
# <http://iopscience.iop.org/article/10.1088/0004-637X/746/1/13/meta>`_,
# `Inglis et al. <http://iopscience.iop.org/article/10.1088/0004-637X/777/1/30/meta>`_)

###############################################################################
# The following image of the flare is now fairly iconic.
aia_cutout03_map = sunpy.map.Map(sample_data.AIA_193_CUTOUT03_IMAGE)
fig = plt.figure()
aia_cutout03_map.plot()
plt.show()

###############################################################################
# Let's take a look at the GOES XRS data.
goes = sunpy.timeseries.TimeSeries(sample_data.GOES_XRS_TIMESERIES)
goes.peek()
plt.show()

###############################################################################
# Next let's investigate the AIA full disk images that are available. Please
# note that these images are not at the full AIA resolution.

aia_131_map = sunpy.map.Map(sample_data.AIA_131_IMAGE)
aia_171_map = sunpy.map.Map(sample_data.AIA_171_IMAGE)
aia_211_map = sunpy.map.Map(sample_data.AIA_211_IMAGE)
aia_335_map = sunpy.map.Map(sample_data.AIA_335_IMAGE)
aia_094_map = sunpy.map.Map(sample_data.AIA_094_IMAGE)
aia_1600_map = sunpy.map.Map(sample_data.AIA_1600_IMAGE)

fig = plt.figure(figsize=(6, 28))
ax = fig.add_subplot(611, projection=aia_131_map)
aia_131_map.plot()
aia_131_map.draw_grid()
ax = fig.add_subplot(612, projection=aia_171_map)
aia_171_map.plot()
aia_171_map.draw_grid()
ax = fig.add_subplot(613, projection=aia_211_map)
aia_211_map.plot()
aia_211_map.draw_grid()
ax = fig.add_subplot(614, projection=aia_335_map)
aia_335_map.plot()
aia_335_map.draw_grid()
ax = fig.add_subplot(615, projection=aia_094_map)
aia_094_map.plot()
aia_094_map.draw_grid()
ax = fig.add_subplot(616, projection=aia_1600_map)
aia_1600_map.plot()
aia_1600_map.draw_grid()
plt.tight_layout(pad=1.50)
plt.show()

###############################################################################
# We also provide a series of AIA cutouts so that you can get a sense of the
# dynamics of the in-falling material.
aia_cutout01_map = sunpy.map.Map(sample_data.AIA_193_CUTOUT01_IMAGE)
aia_cutout02_map = sunpy.map.Map(sample_data.AIA_193_CUTOUT02_IMAGE)
aia_cutout03_map = sunpy.map.Map(sample_data.AIA_193_CUTOUT03_IMAGE)
aia_cutout04_map = sunpy.map.Map(sample_data.AIA_193_CUTOUT04_IMAGE)
aia_cutout05_map = sunpy.map.Map(sample_data.AIA_193_CUTOUT05_IMAGE)

fig = plt.figure(figsize=(6, 28))
ax = fig.add_subplot(511, projection=aia_cutout01_map)
aia_cutout01_map.plot()
ax = fig.add_subplot(512, projection=aia_cutout02_map)
aia_cutout02_map.plot()
ax = fig.add_subplot(513, projection=aia_cutout03_map)
aia_cutout03_map.plot()
ax = fig.add_subplot(514, projection=aia_cutout04_map)
aia_cutout04_map.plot()
ax = fig.add_subplot(515, projection=aia_cutout05_map)
aia_cutout05_map.plot()
plt.tight_layout(pad=2.50)
plt.show()

###############################################################################
# There are a number of other data sources, such as SWAP.
swap_map = sunpy.map.Map(sample_data.SWAP_LEVEL1_IMAGE)
fig = plt.figure()
ax = plt.subplot(projection=swap_map)
swap_map.plot()
plt.show()

###############################################################################
# As well as RHESSI.
rhessi_map = sunpy.map.Map(sample_data.RHESSI_IMAGE)
fig = plt.figure()
#ax = plt.subplot(projection=rhessi_map)
ax = plt.subplot()
rhessi_map.plot()
plt.show()
