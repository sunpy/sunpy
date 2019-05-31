"""
========================
Find peaks in TimeSeries
========================

How to find minimum or maximum peaks in a TimeSeries.
Note: Peak finding is a complex problem that has many potential solutions and
this example is just one method of many.
"""
import numpy as np
import matplotlib.pyplot as plt

from sunpy.timeseries import TimeSeries
from sunpy.data.sample import NOAAINDICES_TIMESERIES as noaa_ind

##############################################################################
# We will now create a TimeSeries object from an observational data source,
# Also, we will truncate it to do analysis on a smaller time duration of 10
# years.

ts_noaa_ind = TimeSeries(noaa_ind, source='NOAAIndices')
my_timeseries = ts_noaa_ind.truncate('1991/01/01', '2001/01/01')
fig, ax = plt.subplots()
my_timeseries.plot()

##############################################################################
# To find extrema in any TimeSeries, we first define a function findpeaks that
# takes in input an iterable data series and a DELTA value. The DELTA value
# controls how much difference between values in the TimeSeries defines an
# extremum point. Inside the function, we iterate over the data values of
# TimeSeries and consider a point to be a local maxima if it has the maximal
# value, and was preceded (to the left) by a value lower by DELTA. Similar
# logic applies to find a local minima.


def findpeaks(series, DELTA):
    """
    Finds extrema in a pandas series data.

    Parameters
    ----------
    series : `pandas.Series`
        The data series from which we need to find extrema.

    DELTA : `float`
        The minimum difference between data values that defines a peak.

    Returns
    -------
    minpeaks, maxpeaks : `list`
        Lists consisting of pos, val pairs for both local minima points and
        local maxima points.
    """
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
minpeaks, maxpeaks = findpeaks(series, DELTA=10.)
# Plotting the figure and extremum points
fig, ax = plt.subplots()
ax.set_ylabel('Sunspot Number')
ax.set_xlabel('Time')
ax.set_title('Peaks in TimeSeries')
series.plot()
ax.scatter(*zip(*minpeaks), color='red', label='min')
ax.scatter(*zip(*maxpeaks), color='green', label='max')
ax.legend()
ax.grid(True)

plt.show()
