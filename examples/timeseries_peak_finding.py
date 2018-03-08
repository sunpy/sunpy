"""
==============================
Find Peaks in SunPy TimeSeries
==============================

This example illustrates how to find minimum or maximum peaks in a TimeSeries.
"""

##############################################################################
# Start by importing the necessary modules.

import numpy as np
import matplotlib.pyplot as plt

import sunpy.data.sample
from sunpy.timeseries import TimeSeries

##############################################################################
# We will now create a TimeSeries object from an observational data source and 
# do analysis on its truncated version.

ts_noaa_ind = sunpy.timeseries.TimeSeries(
    sunpy.data.sample.NOAAINDICES_TIMESERIES, source='NOAAIndices')
my_timeseries = ts_noaa_ind.truncate(0,100)
my_timeseries.peek()

##############################################################################
# Now we take the column 'sunspot SWO' of this TimeSeries and try to find it's
# peaks. For peak finding, we first decide a DELTA value which controls how
# much difference between values in the TimeSeries defines a peak. Then we
# iterate over the data values and consider a point to be a maximum peak if it
# has the maximal value, and was preceded (to the left) by a value lower by 
# DELTA. Similar logic applies to find the minimum valued peaks.

series = my_timeseries.data['sunspot SWO']
DELTA = 10
# Set inital values
mn, mx = np.Inf, -np.Inf
minpos, maxpos = np.NaN, np.NaN
minpeaks = []
maxpeaks = []
lookformax = True
# Iterate over data values
for i in np.arange(len(series)):
    this = series[i]
    if this > mx:
        mx = this
        mxpos = series.index[i]
    if this < mn:
        mn = this
        mnpos = series.index[i]
    if lookformax:
        if this < mx-DELTA:
        	# a local maxima peak
            maxpeaks.append((mxpos, mx))
            mn = this
            mnpos = series.index[i]
            lookformax = False
    else:
        if this > mn+DELTA:
        	# a local minima peak
            minpeaks.append((mnpos, mn))
            mx = this
            mxpos = series.index[i]
            lookformax = True
# Plotting the figure and peaks
plt.figure(2)
fontsize = 14
plt.ylabel(r'Sunspot Number')
plt.xlabel(r'Time')
plt.title(r'Peaks in TimeSeries')
plt.tight_layout()
series.plot(legend=True)
plt.scatter(np.array(minpeaks)[:,0],np.array(minpeaks)[:,1],color='red', label='min')
plt.scatter(np.array(maxpeaks)[:,0],np.array(maxpeaks)[:,1],color='green', label='max')
plt.legend()
plt.grid(True)
plt.show()