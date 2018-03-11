"""
==============================
Find Peaks in SunPy TimeSeries
==============================

This example illustrates how to find minimum or maximum peaks in a TimeSeries.
Note: Peak finding is a complex problem that has many potential solutions and
this example is just one method of many.
"""

##############################################################################
# Start by importing the necessary modules.

import numpy as np
import matplotlib.pyplot as plt

import sunpy.data.sample
from sunpy.timeseries import TimeSeries

##############################################################################
# We will now create a TimeSeries object from an observational data source,
# Also, we will truncate it to do analysis on a smaller time duration of 10
# years.

ts_noaa_ind = sunpy.timeseries.TimeSeries(
    sunpy.data.sample.NOAAINDICES_TIMESERIES, source='NOAAIndices')
my_timeseries = ts_noaa_ind.truncate('1991/01/01', '2001/01/01')
my_timeseries.peek()

##############################################################################
# Now we take the column 'sunspot SWO' of this TimeSeries and try to find it's
# extrema. To do this, we first decide a DELTA value which controls how much
# difference between values in the TimeSeries defines an extremum point. Then
# we iterate over the data values and consider a point to be a local maxima if
# it has the maximal value, and was preceded (to the left) by a value lower by
# DELTA. Similar logic applies to find a local minima.

series = my_timeseries.data['sunspot SWO']
# We take delta to be approximately the length of smallest peak that we wish
# to detect
DELTA = 10
# Set inital values
mn, mx = np.Inf, -np.Inf
minpeaks = []
maxpeaks = []
lookformax = True
# Iterate over items in series
for time_pos, value in series.iteritems():
    this = value
    if this > mx:
        mx = this
        mxpos = time_pos
    if this < mn:
        mn = this
        mnpos = time_pos
    if lookformax:
        if this < mx-DELTA:
            # a local maxima
            maxpeaks.append((mxpos, mx))
            mn = this
            mnpos = time_pos
            lookformax = False
    else:
        if this > mn+DELTA:
            # a local minima
            minpeaks.append((mnpos, mn))
            mx = this
            mxpos = time_pos
            lookformax = True
# Plotting the figure and extremum points
plt.figure()
plt.ylabel('Sunspot Number')
plt.xlabel('Time')
plt.title('Peaks in TimeSeries')
series.plot()
plt.scatter(*zip(*minpeaks), color='red', label='min')
plt.scatter(*zip(*maxpeaks), color='green', label='max')
plt.legend()
plt.grid(True)
plt.show()
