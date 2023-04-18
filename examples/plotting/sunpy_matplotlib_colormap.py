"""
=========================================
Using the sunpy Colormaps with matplotlib
=========================================

How you can use the sunpy colormaps with matplotlib.
"""
import matplotlib
import matplotlib.pyplot as plt
import numpy as np

import sunpy.visualization.colormaps as cm

###############################################################################
# When the sunpy colormaps are imported, the sunpy colormaps are registered
# with matplotlib. It is now possible to access the colormaps with the following command.

sdoaia171 = matplotlib.colormaps['sdoaia171']

###############################################################################
# You can get the list of all sunpy colormaps with:

print(cm.cmlist.keys())

###############################################################################
# Let's now create a data array.

delta = 0.025
x = y = np.arange(-3.0, 3.0, delta)
X, Y = np.meshgrid(x, y)
Z1 = np.exp(-X**2 - Y**2)
Z2 = np.exp(-(X - 1)**2 - (Y - 1)**2)
Z = (Z1 - Z2) * 2

###############################################################################
# Let's now plot the results with the colormap.

fig, ax = plt.subplots()
im = ax.imshow(
    Z, interpolation='bilinear', cmap=sdoaia171,
    origin='lower', extent=[-3, 3, -3, 3],
    vmax=abs(Z).max(), vmin=-abs(Z).max()
)
plt.show()
