"""
======================
Masking the solar disk
======================

This example shows how to mask off emission from the disk.
"""
import numpy as np
import matplotlib.pyplot as plt
import numpy.ma as ma
from sunpy.data.sample import AIA_171_IMAGE
from sunpy.map import Map
###############################################################################
# We first create the Map using the sample data.
aia = Map(AIA_171_IMAGE)

###############################################################################
# Next we build two arrays which include all of the x and y pixel indices.
# We must not forget to add the correct units because we will next pass this_sji_map
# into a sunpy function which requires them.
x, y = np.meshgrid(*[np.arange(v.value) for v in aia.dimensions]) * u.pixel

###############################################################################
# Now we can convert this to helioprojective coordinates and create a new
# array which contains the normalized radial position for each pixel
hpc_x, hpc_y = aia.pixel_to_data(x,y)
r = np.sqrt(hpc_x ** 2 + hpc_y ** 2) / aia.rsun_obs

###############################################################################
# Finally, we create a mask where all values which are less then Rsun are
# masked. We also make a slight change to the colormap so that masked values
# are shown as black (instead of the default white).
mask = ma.masked_less_equal(r, 1)
palette = aia.plot_settings['cmap']
palette.set_bad('black')

###############################################################################
# Now we create a new custom aia with our new mask and
# plot the result using our modified colormap
scaled_map = Map(aia.data, aia.meta, mask=mask.mask)
scaled_map.plot(cmap=palette)
scaled_map.draw_limb()
plt.show()
