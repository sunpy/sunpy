"""
==================================================
Example showcase obtaining and plotting METIS data
==================================================

This example demonstrates how to download, load, and visualize visible-light
total brightness observations from the METIS coronagraph aboard
the Solar Orbiter mission.

METIS observes the solar corona by occulting the bright solar disk,
allowing faint coronal structures to be imaged in visible light and
ultraviolet wavelengths.
"""
# sphinx_gallery_tags = ["Map", "METIS", "SOAR", "Solar Orbiter", "Visualization"]
# sphinx_gallery_thumbnail_number = -1

import matplotlib.pyplot as plt
import numpy as np
from matplotlib import colors

import sunpy.map
from sunpy.net import Fido
from sunpy.net import attrs as a

###############################################################################
# First, download some METIS data with `~sunpy.net.Fido`.

result = Fido.search(a.Time("2022-03-22 21:00", "2022-03-22 22:50"),
                     a.Instrument.metis, a.Level(2), a.Provider.soar,
                     a.soar.Product.metis_vl_tb)
metis_data_product = Fido.fetch(result)

###############################################################################
# The downloaded METIS FITS file contains multiple image header data units (HDUs).
# When passed to `sunpy.map.Map`, SunPy creates a MapSequence containing one map
# for each HDU.
# We then filter the sequence to obtain only the VL total-brightness maps.
# Finally, we select the first VL total-brightness map.

vl_maps = sunpy.map.Map(metis_data_product[0])
vl_tb_list = [m for m in vl_maps if m.measurement == "VL-TB"]
vl_tb_example = vl_tb_list[0]

###############################################################################
# A `~sunpy.map.sources.METISMap` is created with a default mask property that
# excludes pixels inside the inner occulter and outside the outer field of view

fig = plt.figure()
ax = fig.add_subplot(projection=vl_tb_example)
im = vl_tb_example.plot(axes=ax)
fig.colorbar(im)

##############################################################################
# In this example, we also know there is an eruption event.
# We start by animating this time series.
# To do this we load all the HDUs of all the files and filter down to
# just the "VL-TB" HDUs and make a `~sunpy.map.MapSequence`.

vl_tb_maps = sunpy.map.Map(metis_data_product)
m_seq_vltb = sunpy.map.Map([m for m in vl_tb_maps if m.measurement == "VL-TB"], sequence=True)
fig = plt.figure()
m_seq_vltb.plot()

###############################################################################
# To get a better idea of what the event looks like we need to create
# a running difference.

m_seq_running = sunpy.map.Map(
    [m - prev_m.quantity for m, prev_m in zip(m_seq_vltb[1:], m_seq_vltb[:-1])],
    sequence=True
)

###############################################################################
# Now we can normalise the values from all the maps
# and create a running difference animation

all_diff_data = np.concatenate([m.data for m in m_seq_running])
vmin, vmax = np.percentile(all_diff_data, [0.05, 99.95])
norm = colors.Normalize(vmin=vmin, vmax=vmax)

fig = plt.figure()
ax = fig.add_subplot(projection=m_seq_running.maps[0])
ani_running = m_seq_running.plot(
    axes=ax,
    title='VL-TB Running Difference',
    norm=norm
)
plt.colorbar(extend='both', label=m_seq_running[0].unit.to_string())
plt.show()
