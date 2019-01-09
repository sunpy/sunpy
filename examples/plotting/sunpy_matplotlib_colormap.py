"""
========================================
Using the SunPy Colormaps and Matplotlib
========================================

This examples shows how you can use the SunPy colormaps with Matplotlib.
Also the full range of colormaps we provide.
"""
###############################################################################
# When the sunpy colormaps are imported with the following command, we register the
# SunPy colormap names with Matplotlib.

import sunpy.cm as cm

###############################################################################
# You can now get at any of the colormaps, for example 'sdoaia171', with the
# following command

import matplotlib.pyplot as plt
sdoaia171 = plt.get_cmap('sdoaia171')

###############################################################################
# We can now create a standard matplotlib plot.

import numpy as np
delta = 0.025
x = y = np.arange(-3.0, 3.0, delta)
X, Y = np.meshgrid(x, y)
Z1 = np.exp(-X**2 - Y**2)
Z2 = np.exp(-(X - 1)**2 - (Y - 1)**2)
Z = (Z1 - Z2) * 2

fig, ax = plt.subplots()
im = ax.imshow(Z, interpolation='bilinear', cmap=sdoaia171,
               origin='lower', extent=[-3, 3, -3, 3],
               vmax=abs(Z).max(), vmin=-abs(Z).max())

plt.show()

###############################################################################
# If you don't remember what colormaps are available, you can get the list with:

print(cm.cmlist.keys())

###############################################################################
# We also provide a function that will display all our colormaps:

# The next line sets the thumbnail for the second figure in the gallery.
# sphinx_gallery_thumbnail_number = 2
cm.show_colormaps()
