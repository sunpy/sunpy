"""
==============================
Plotting the solar cycle index
==============================

This example demonstrates how to plot a solar cycle in terms of
the number of sunspots and a prediction for the next few years.
"""
import datetime
import matplotlib.pyplot as plt

import sunpy.timeseries as ts
from sunpy.time import parse_time
from sunpy.net import Fido, attrs as a

###############################################################################
# For this example we will downaload sunpy's data which has an old
# sunspot record from 2015 and prediction for the cycle till 2019 from (SWPC)
# part of the National Oceanic and Atmospheric Administration (NOAA).
# This code snippet downloads and grabs the NOAA solar cycle data as a
# ``TimeSeries``.
result = Fido.search(a.Time("2015/1/1", datetime.datetime.now()),
                     a.Instrument('noaa-indices'))
f_noaa_indices = Fido.fetch(result)
result = Fido.search(a.Time(datetime.datetime.now(), datetime.datetime.now() +
                     datetime.timedelta(days=4*365)),
                     a.Instrument(''noaa-predict''))
f_noaa_predict = Fido.fetch(result)

noaa = ts.TimeSeries(f_noaa_indices, source='noaaindices')
noaa_predict = ts.TimeSeries(f_noaa_predict, source='noaapredictindices')

###############################################################################
# Finally, we plot both ``noaa`` and ``noaa_predict`` together, with an arbitrary
# range for the strength of the next solar cycle.
# Here ``noaa`` refers to existing solar cycle data from 2015 and
# ``noaa_predict`` is a prediction for the rest of the solar cycle till 2019.
fig, ax = plt.subplots()
ax.plot(noaa.data.index, noaa.data['sunspot RI'], label='Sunspot Number')
ax.plot(noaa_predict.data.index, noaa_predict.data['sunspot'],
        color='grey', label='Near-term Prediction')
ax.fill_between(noaa_predict.data.index, noaa_predict.data['sunspot low'],
                noaa_predict.data['sunspot high'], alpha=0.3, color='grey')
ax.set_ylim(bottom=0)
ax.set_ylabel('Sunspot Number')
ax.set_xlabel('Year')
ax.legend()

plt.show()
