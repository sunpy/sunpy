# coding: utf-8
"""
===============================================
Editing the colormap and normalization of a Map
===============================================

How to edit the display of a map.
"""
import matplotlib.pyplot as plt
import matplotlib.colors as colors

import sunpy.map
from sunpy.data.sample import AIA_171_IMAGE

###############################################################################
# We start with the sample data
aiamap = sunpy.map.Map(AIA_171_IMAGE)

###############################################################################
# How a Map is displayed is determined by its colormap, which sets the colors
# , and the normalization, which sets how data values are translated to colors.
# Lets set the colormap and normalization, to replace them later.
cmap = plt.get_cmap('Greys_r')
norm = colors.LogNorm(100, aiamap.max())

###############################################################################
# To see all of the colormaps SunPy provides see `sunpy.visualization.colormaps`.
# Matplotlib provides a number of `colormaps <https://matplotlib.org/examples/color/colormaps_reference.html>`_
# and `normalizations <https://matplotlib.org/users/colormapnorms.html>`_.
# For more advanced normalizations see `astropy.visualization`.
# The colormap and normalization for a map are passed on as attributes to `plot`
ax = plt.subplot(projection=aiamap)
aiamap.plot(cmap=cmap, norm=norm)
plt.colorbar()
plt.show()
