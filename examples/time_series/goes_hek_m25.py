"""
==============================
Flare times on a GOES XRS plot
==============================

How to plot flare times as provided by the HEK on a GOES XRS plot.
"""
import matplotlib.pyplot as plt

from sunpy.net import Fido
from sunpy.net import attrs as a
from sunpy.time import parse_time
from sunpy.timeseries import TimeSeries

###############################################################################
# Let's grab GOES XRS data for a particular time of interest and the HEK flare
# data for this time from the NOAA Space Weather Prediction Center (SWPC)

tr = a.Time('2011-06-07 04:00', '2011-06-07 12:00')
results = Fido.search(tr, a.Instrument.xrs | a.hek.FL & (a.hek.FRM.Name == 'SWPC'))

###############################################################################
# Then download the XRS data and load it into a TimeSeries
files = Fido.fetch(results)
goes = TimeSeries(files)

###############################################################################
# Next lets load the HEK metadata from Fido result to a HEKRow
hek_results = results.get_response(1)
flares_hek = hek_results[0]

###############################################################################
# Lets plot everything together
fig, ax = plt.subplots()
goes.plot()
ax.axvline(parse_time(flares_hek['event_peaktime']).datetime)
ax.axvspan(parse_time(flares_hek['event_starttime']).datetime,
           parse_time(flares_hek['event_endtime']).datetime,
           alpha=0.2, label=flares_hek['fl_goescls'])
ax.legend(loc=2)
ax.set_yscale('log')
plt.show()
