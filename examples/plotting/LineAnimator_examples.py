"""
===========================
How to use the LineAnimator
===========================

This example shows off some ways in which you can use the
LineAnimator object to animate line plots.
"""
import matplotlib.pyplot as plt
import numpy as np

from sunpy.visualization.animator import LineAnimator

###############################################################################
# Animate a 2D cube of random data as a line plot along an
# axis where the x-axis drifts with time.

# Define some random data
data_shape0 = (10, 20)
data0 = np.random.rand(*data_shape0)

###############################################################################
# Define the axis that will make up the line plot
plot_axis0 = 1
slider_axis0 = 0

###############################################################################
# Let's customize the values along the x-axis.  To do this, we must define the
# edges of the pixels/bins being plotted along the x-axis.  This requires us to
# supply an array, say xdata, of length equal to data.shape[plot_axis_index]+1.
# In this example, the data has a shape of (10, 20) and let's say we are
# iterating through the 0th axis and plotting the 1st axis,
# i.e. plot_axis_index=1.  Therefore we need to define an xdata array of length
# 21.
# This will give the same customized x-axis values for each frame of the
# animation. However, what if we want the x-axis values to change as we
# animate through the other dimensions of the cube?  To do this we supply a
# (10, 21) xdata where each row (i.e. xdata[i, :]) gives the pixel/bin edges
# along the x-axis for the of the i-th frame of the animation.  Note that this
# API extends in the same way to higher dimension.  In our 2D case here though,
# we can define our non-constant x-axis values like so:
xdata = np.tile(np.linspace(0, 100, (data_shape0[plot_axis0] + 1)), (data_shape0[slider_axis0], 1))

###############################################################################
# Generate animation object with variable x-axis data.
ani = LineAnimator(data0, plot_axis_index=plot_axis0, axis_ranges=[None, xdata])

###############################################################################
# Show plot
plt.show()
