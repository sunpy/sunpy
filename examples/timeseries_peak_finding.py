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
# To find extrema in any TimeSeries, we first define a function findpeaks that
# takes in input an iterable data series and a DELTA value. The DELTA value
# controls how much difference between values in the TimeSeries defines an
# extremum point. Inside the function, we iterate over the data values of
# TimeSeries and consider a point to be a local maxima if it has the maximal
# value, and was preceded (to the left) by a value lower by DELTA. Similar
# logic applies to find a local minima.


def findpeaks(series, DELTA):
    # Set inital values
    mn, mx = np.Inf, -np.Inf
    minpeaks = []
    maxpeaks = []
    lookformax = True
    start = True
    # Iterate over items in series
    for time_pos, value in series.iteritems():
        if value > mx:
            mx = value
            mxpos = time_pos
        if value < mn:
            mn = value
            mnpos = time_pos
        if lookformax:
            if value < mx-DELTA:
                # a local maxima
                maxpeaks.append((mxpos, mx))
                mn = value
                mnpos = time_pos
                lookformax = False
            elif start:
                # a local minima at beginning
                minpeaks.append((mnpos, mn))
                mx = value
                mxpos = time_pos
                start = False
        else:
            if value > mn+DELTA:
                # a local minima
                minpeaks.append((mnpos, mn))
                mx = value
                mxpos = time_pos
                lookformax = True
    # check for extrema at end
    if value > mn+DELTA:
        maxpeaks.append((mxpos, mx))
    elif value < mx-DELTA:
        minpeaks.append((mnpos, mn))
    return minpeaks, maxpeaks

##############################################################################
# Now we take the column 'sunspot SWO' of this TimeSeries and try to find it's
# extrema using the function findpeaks. We take the value of DELTA to be
# approximately the length of smallest peak that we wish to detect.

series = my_timeseries.data['sunspot SWO']
minpeaks, maxpeaks = findpeaks(series, DELTA=10)
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
