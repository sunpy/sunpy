"""
=======================================
Editing Map Colormaps and Normalization
=======================================

A simple example to show how to edit the display of a map
"""

import matplotlib.cm as cm
import matplotlib.colors as colors
import matplotlib.pyplot as plt
import sunpy.map
from sunpy.data.sample import AIA_171_IMAGE


###############################################################################
# We first create the Map using the sample data.

aiamap = sunpy.map.Map(AIA_171_IMAGE)

###############################################################################
# Now lets replace the colormap which sets the colors as well as the
# normalization which sets how data values are translated to colors

aiamap.plot_settings['cmap'] = cm.get_cmap('Greys_r')
aiamap.plot_settings['norm'] = colors.LogNorm(100, aiamap.max())

###############################################################################
# You can find more colormaps in matplotlib (http://matplotlib.org/examples/color/colormaps_reference.html)
# or look at the sunpy colormaps in `sunpy.cm`
# For more normalizations check out `matplotlib.colors.Normalize`
# or astropy provides additional norms in `astropy.visualization`

aiamap.plot()
plt.colorbar()
plt.show()
