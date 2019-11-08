"""
Fading between two maps
=======================

How to plot two maps on top of each other, with a slider to fade between them.
"""

import matplotlib.pyplot as plt
from matplotlib.widgets import Slider

import sunpy.map
from sunpy.data.sample import AIA_171_IMAGE, AIA_1600_IMAGE

###############################################################################
# Start by loading two AIA maps from the sample data
map_171 = sunpy.map.Map(AIA_171_IMAGE)
map_1600 = sunpy.map.Map(AIA_1600_IMAGE)

###############################################################################
# Create a figure, and add the axes that will show the maps
fig = plt.figure()
ax = fig.add_axes([0.1, 0.2, 0.9, 0.8], projection=map_171)

###############################################################################
# Plot both maps on the same axes, and store the image objects
im_1600 = map_1600.plot(axes=ax)
im_171 = map_171.plot(axes=ax, alpha=0.5)
ax.set_title('AIA 171 + 1600')

###############################################################################
# Now add another axes, add add a slider to it
ax_slider = fig.add_axes([0.25, 0.1, 0.65, 0.03])
slider = Slider(ax_slider, 'Alpha', 0, 1, valinit=0.5)


###############################################################################
# Finally, define what happens when the slider is changed and link this to the
# slider
def update(val):
    alpha = slider.val
    im_171.set_alpha(alpha)


slider.on_changed(update)

plt.show()
