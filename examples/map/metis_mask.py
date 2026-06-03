"""
==============
Plotting METIS
==============

In this example, we plot the VL total brightness data product from METIS with the recommended contrast cutoffs.
"""
# sphinx_gallery_tags = ["Map", "METIS"]

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
# Then we use one of the data products to create a simple map.
# A `~sunpy.map.sources.METISMap` is created with a default mask property that
# excludes pixels inside the inner occulter and outside the outer field of view

vl_maps = sunpy.map.Map(metis_data_product[0])
vl_tb_list = [m for m in vl_maps if m.measurement == "VL-TB"]
vl_tb_example = vl_tb_list[0]

fig = plt.figure()
ax = fig.add_subplot(projection=vl_tb_example)
norm = vl_tb_example.plot_settings["norm"]
vl_tb_example.plot(axes=ax, cmap=vl_tb_example.plot_settings["cmap"])
plt.show()



##############################################################################
# In this example, we also know there is a related event. To get a better idea
# of what the event looks like we need to create a running difference.
# To do this we need to create a list of maps and filter for the VL-TB.


vl_tb_maps = sunpy.map.Map(metis_data_product,sequence=True)
m_seq_vltb = [m for m in vl_tb_maps if m.measurement == "VL-TB"]


###############################################################################
# Now we calculate the running difference.

m_seq_running = sunpy.map.Map(
    [m - prev_m.quantity for m, prev_m in zip(m_seq_vltb[1:], m_seq_vltb[:-1])],
    sequence=True
)


###############################################################################
# Now we can normalise the values from all the maps
# and create a running difference animation

all_diff_data = np.concatenate([m.data for m in m_seq_running])
vmin, vmax = np.percentile(all_diff_data, [1, 99])
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
