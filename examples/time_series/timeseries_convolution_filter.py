"""
======================================================
Smoothing of timeSeries data using convolution filters
======================================================

How to smooth a TimeSeries using a convolution filter
kernel from `~astropy.convolution` and `~astropy.convolution.convolve`
function.
"""
import matplotlib.pyplot as plt

from astropy.convolution import Box1DKernel, convolve

from sunpy.data.sample import NOAAINDICES_TIMESERIES as noaa_ind
from sunpy.timeseries import TimeSeries

###############################################################################
# Let's first create a TimeSeries from sample data
ts_noaa_ind = TimeSeries(noaa_ind, source='NOAAIndices')

###############################################################################
# Now we will extract data values from the TimeSeries and apply a BoxCar filter
# to get smooth data. Boxcar smoothing is equivalent to taking our signal and
# using it to make a new signal where each element is the average of w adjacent
# elements. Here we will use AstroPy’s convolve function with a “boxcar” kernel
# of width w = 10.
ts_noaa_ind = ts_noaa_ind.add_column(
    'sunspot SWO Smoothed',
    convolve(ts_noaa_ind.quantity('sunspot SWO'), kernel=Box1DKernel(10))
)

###############################################################################
# Plotting original and smoothed timeseries
plt.ylabel('Sunspot Number')
plt.xlabel('Time')
plt.title('Smoothing of Time Series')
plt.plot(ts_noaa_ind.quantity('sunspot SWO'), label='original data')
plt.plot(ts_noaa_ind.quantity('sunspot SWO Smoothed'), label='smoothed')
plt.legend()
plt.show()
