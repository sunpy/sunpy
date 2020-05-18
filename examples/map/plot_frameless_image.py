"""
===============================
Plotting a Map without any Axes
===============================

This examples shows you how to plot a Map without any annotations at all, i.e.
to save as an image.
"""
import matplotlib.pyplot as plt

##############################################################################
# Start by importing the necessary modules.
import astropy.units as u

import sunpy.map
from sunpy.data.sample import AIA_171_IMAGE

##############################################################################
# Create a `sunpy.map.GenericMap`.
smap = sunpy.map.Map(AIA_171_IMAGE)

##############################################################################
# Plot the Map without a frame.

# Setup a frameless figure and an axes which spans the whole canvas.
figure = plt.figure(frameon=False)
axes = plt.Axes(figure, [0., 0., 1., 1.])

# Disable the axis and add them to the figure.
axes.set_axis_off()
figure.add_axes(axes)

# Plot the map without any annotations
# This might raise a warning about the axes being wrong but we can ignore this
# as we are not plotting any axes.
im = smap.plot(axes=axes, annotate=False, clip_interval=(1, 99.99)*u.percent)

##############################################################################
# At this point you could save the figure with ``plt.savefig()`` or show it:

plt.show()
