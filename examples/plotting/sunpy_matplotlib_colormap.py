"""
========================================
Using the SunPy Colormaps and matplotlib
========================================

This examples shows how you can use the SunPy colormaps with Matplotlib.
"""
###############################################################################
# When import the sunpy colormaps with the following command we register the
# SunPy colormap names with mapltotlib.

import matplotlib.pyplot as plt
import sunpy.cm

###############################################################################
# You can now get at any of the colormaps, for example 'sdoaia171', with the
# following command

sdoaia171 = plt.get_cmap('sdoaia171')

###############################################################################
# Let's now make a standard matplotlib plot

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
# If you don't remember what colormaps are available, you can get the list with
print(sunpy.cm.cmlist.keys())
