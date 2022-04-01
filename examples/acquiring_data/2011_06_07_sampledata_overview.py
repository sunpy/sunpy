"""
========================
Sample data set overview
========================

An overview of the coordinated sample data set (available through `sunpy.data.sample`).
"""
import matplotlib.pyplot as plt
import numpy as np

import astropy.units as u
from astropy.coordinates import SkyCoord

import sunpy.coordinates
import sunpy.data.sample as sample_data
import sunpy.map
import sunpy.timeseries
from sunpy.io.special import srs

###############################################################################
# On 2011 June 7, various solar instruments observed a spectacular solar
# eruption from NOAA AR 11226. The event included an M2.5 flare, a
# filament eruption, a coronal mass ejection, and a global coronal EUV wave (IAU standard:
# SOL2011-06-07T06:24:00L045C112). This event was spectacular because it
# features the ejection of a large amount of prominence material, much of which
# failed to escape and fell back to the solar surface.
# This event received some press coverage (e.g. `National Geographics
# <https://news.nationalgeographic.com/news/2011/06/110608-solar-flare-sun-science-space/>`_,
# `Discover Magazine <http://blogs.discovermagazine.com/badastronomy/2011/06/07/the-sun-lets-loose-a-huge-explosion/>`_)
# and the literature contains a number of a papers about it (e.g. `Li et al.
# <https://iopscience.iop.org/article/10.1088/0004-637X/746/1/13/meta>`_,
# `Inglis et al. <https://iopscience.iop.org/article/10.1088/0004-637X/777/1/30/meta>`_)
# The following image of the flare is now fairly iconic.

aia_cutout03_map = sunpy.map.Map(sample_data.AIA_193_CUTOUT03_IMAGE)
plt.figure()
aia_cutout03_map.plot()
plt.show()

###############################################################################
# Let's take a look at the GOES XRS data.

goes = sunpy.timeseries.TimeSeries(sample_data.GOES_XRS_TIMESERIES)
fig = plt.figure()
goes.plot()
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
for i, m in enumerate([aia_131_map, aia_171_map, aia_211_map,
                       aia_335_map, aia_094_map, aia_1600_map]):
    ax = fig.add_subplot(6, 1, i+1, projection=m)
    m.plot(axes=ax, clip_interval=(0.5, 99.9)*u.percent)
    m.draw_grid(axes=ax)
fig.tight_layout(pad=8.50)
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
for i, m in enumerate([aia_cutout01_map, aia_cutout02_map, aia_cutout03_map,
                       aia_cutout04_map, aia_cutout05_map]):
    ax = fig.add_subplot(5, 1, i+1, projection=m)
    m.plot(axes=ax)
fig.tight_layout(pad=5.50)
plt.show()

###############################################################################
# There are a number of other data sources available as well, such as SWAP.

swap_map = sunpy.map.Map(sample_data.SWAP_LEVEL1_IMAGE)
plt.figure()
swap_map.plot()
plt.show()

###############################################################################
# Also RHESSI.

rhessi_map = sunpy.map.Map(sample_data.RHESSI_IMAGE)
plt.figure()
rhessi_map.plot()
plt.show()

##############################################################################
# NORH

norh = sunpy.timeseries.TimeSeries(sample_data.NORH_TIMESERIES)
fig = plt.figure()
norh.plot()
plt.show()

####################################################################################
# NOAA overlaid with HMI

noaa = srs.read_srs(sample_data.SRS_TABLE)
smap = sunpy.map.Map(sample_data.HMI_LOS_IMAGE)

noaa = noaa[np.logical_or(noaa['ID'] == 'I', noaa['ID'] == 'IA')]

lats = noaa['Latitude']
lngs = noaa['Longitude']
numbers = noaa['Number']

fig = plt.figure()
ax = fig.add_subplot(projection=smap)

smap.plot_settings["norm"].vmin = -200
smap.plot_settings["norm"].vmax = 200
smap.plot(axes=ax)
smap.draw_limb(axes=ax)

c = SkyCoord(lngs, lats, frame="heliographic_stonyhurst")

ax.plot_coord(c, 'o')

for lat, lng, num in zip(lats.value, lngs.value, numbers):
    ax.annotate(num, (lng, lat),
                xycoords=ax.get_transform('heliographic_stonyhurst'),
                color='blue',
                fontweight='bold')

plt.show()

####################################################################################
# EVE

eve = sunpy.timeseries.TimeSeries(sample_data.EVE_TIMESERIES, source='EVE')
fig = plt.figure(figsize=(10, 6))
eve.plot()
plt.legend(bbox_to_anchor=(1.01, 0.85), loc='upper left', borderaxespad=0)
fig.tight_layout()
plt.show()

####################################################################################
# LYRA

lyra = sunpy.timeseries.TimeSeries(sample_data.LYRA_LEVEL3_TIMESERIES, source='lyra')
plt.rcParams['figure.figsize'] = [12, 7]
lyra.plot()
plt.show()

####################################################################################
# GBM

gbm = sunpy.timeseries.TimeSeries(sample_data.GBM_TIMESERIES, source='GBMSummary')
fig = plt.figure(figsize=(12, 6))
gbm.plot()
plt.legend(bbox_to_anchor=(1.01, 0.85), loc='upper left', borderaxespad=0)
fig.tight_layout()
plt.show()
