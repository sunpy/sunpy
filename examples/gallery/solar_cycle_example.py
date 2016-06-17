"""
===============
The Solar Cycle
===============

This example shows the current and possible next solar cycle.
"""
import datetime
import matplotlib.pyplot as plt

import sunpy.lightcurve as lc
from sunpy.data.sample import NOAAINDICES_LIGHTCURVE, NOAAPREDICT_LIGHTCURVE

###############################################################################
# For this example we will use the SunPy sample data, if you want the current
# data, delete the argument to the ``create`` function. i.e.
# ``noaa = lc.NOAAIndicesLightCurve.create()``
noaa = lc.NOAAIndicesLightCurve.create(NOAAINDICES_LIGHTCURVE)
noaa_predict = lc.NOAAPredictIndicesLightCurve.create(NOAAPREDICT_LIGHTCURVE)

###############################################################################
# Next lets grab the data again to create a new data structure that we will
# shift by 12 years to simulate the next solar cycle. We will truncate the
# data to only plot what is necessary.
noaa2 = lc.NOAAIndicesLightCurve.create(NOAAINDICES_LIGHTCURVE)
noaa2.data = noaa2.data.shift(2, freq=datetime.timedelta(days=365*12))
noaa2 = noaa2.truncate('2021/04/01', '2030/01/01')

###############################################################################
# Finally lets plot everything together with some arbitrary range for the
# strength of the next solar cycle.
plt.plot(noaa.data.index, noaa.data['sunspot RI'], label='Sunspot Number')
plt.plot(noaa_predict.data.index, noaa_predict.data['sunspot'],
         color='grey', label='Near-term Prediction')
plt.fill_between(noaa_predict.data.index, noaa_predict.data['sunspot low'],
                 noaa_predict.data['sunspot high'], alpha=0.3, color='grey')

plt.fill_between(noaa2.data.index, noaa2.data['sunspot RI smooth']*0.4,
                 noaa2.data['sunspot RI smooth']*1.3, alpha=0.3, color='grey',
                 label='Next Cycle Predict')
plt.ylim(0)
plt.text('2011-01-01', 120, 'Cycle 24', fontsize=16)
plt.text('2024-01-01', 120, 'Cycle 25', fontsize=16)
plt.ylabel('Sunspot Number')
plt.xlabel('Year')
plt.legend(loc=2, framealpha=0.5)
plt.show()
