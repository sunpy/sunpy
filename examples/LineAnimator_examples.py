"""
=============
LineAnimator
=============

This example shows off some ways in which you can use the
LineAnimator object to animate line plots.
"""
import numpy as np
import matplotlib.pyplot as plt

from sunpy.visualization.imageanimator import LineAnimator

# Example 1: Animate a 2D cube of random data as a line plot along an
# axis where the x-axis drifts with time.

# Define some random data
data_shape0 = (10, 20)
data0 = np.random.rand(*data_shape0)

# Define the axis that will make up the line plot
plot_axis0 = 1
slider_axis0 = 0

# Define value along x axis which drift with time.  To do this, define
# xdata to be the same shape as the data where each row/column
# (depending on axis to be animated) represents the x-axis values for
# a single frame of the animations.
xdata = np.tile(np.linspace(0, 100, data_shape0[plot_axis0]), (data_shape0[slider_axis0], 1))

# Generate animation object with variable x-axis data.
ani = LineAnimator(data0, plot_axis_index=plot_axis0, axis_ranges=[None, xdata])

# Show plot
plt.show()
