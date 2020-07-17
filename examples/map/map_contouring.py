# coding: utf-8
"""
=========================
Finding contours of a map
=========================

This example shows how to find and plot contours on a map.
"""

import matplotlib.pyplot as plt

import sunpy.map
from sunpy.data.sample import AIA_193_IMAGE

###############################################################################
# Start by loading the sample data
aiamap = sunpy.map.Map(AIA_193_IMAGE)

###############################################################################
# Now find the contours
contours = aiamap.contour(50000)

##############################################################################
# Finally, plot the map and add the contours
plt.figure()
ax = plt.subplot(projection=aiamap)
aiamap.plot()
for contour in contours:
    ax.plot_coord(contour)
plt.show()
