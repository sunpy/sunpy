"""
==============================
Plotting the solar cycle index
==============================

How to plot the current and possible next solar cycle.
"""
import matplotlib.pyplot as plt

import sunpy.timeseries as ts
from sunpy.data.sample import NOAAINDICES_TIMESERIES, NOAAPREDICT_TIMESERIES

###############################################################################
# For this example we will use the SunPy sample data. This code snippet grabs
# the most current NOAA solar cycle data as a ``TimeSeries``.
noaa = ts.TimeSeries(NOAAINDICES_TIMESERIES, source='noaaindices')
noaa_predict = ts.TimeSeries(NOAAPREDICT_TIMESERIES, source='noaapredictindices')

###############################################################################
# Finally, we plot both ``noaa`` and ``noaa2`` together, with an arbitrary
# range for the strength of the next solar cycle.
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
