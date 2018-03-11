"""
======================================================
Smoothing of TimeSeries Data Using Convolution Filters
======================================================

This example illustrates smoothing a TimeSeries using a convolution filter
kernel from `~astropy.convolution` and `~astropy.convolution.convolve`
function.
"""

##############################################################################
# Start by importing the necessary modules.

import datetime

import numpy as np
import matplotlib.pyplot as plt
from pandas import DataFrame
from astropy.convolution import convolve, Box1DKernel

from sunpy.timeseries import TimeSeries
import sunpy.data.sample

###############################################################################
# Let's first create a TimeSeries from sample data
ts_noaa_ind = sunpy.timeseries.TimeSeries(
    sunpy.data.sample.NOAAINDICES_TIMESERIES, source='NOAAIndices')

# #############################################################################
# Now we will extract data values from the TimeSeries and apply a BoxCar filter
# to get smooth data. Boxcar smoothing is equivalent to taking our signal and
# using it to make a new signal where each element is the average of w adjacent
# elements. Here we will use AstroPy’s convolve function with a “boxcar” kernel
# of width w = 10.

data = ts_noaa_ind.data['sunspot SWO'].values
index = ts_noaa_ind.data['sunspot SWO'].index
# Apply convolution filter
convolved_data = convolve(data, kernel=Box1DKernel(10))
# Plotting original and smoothed timeseries
plt.ylabel('Sunspot Number')
plt.xlabel('Time')
plt.title('Smoothing of Time Series')
plt.plot(index, data, label='Sunspot SWO')
plt.plot(index, convolved_data, label='Sunspot SWO Smoothed')
plt.legend()
plt.show()
