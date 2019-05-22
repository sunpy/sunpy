"""
===============================
GOES Flare and HEK Plot Example
===============================

An example showing how to combine GOES and HEK data
"""
import matplotlib.pyplot as plt

from sunpy.timeseries import TimeSeries
from sunpy.time import TimeRange, parse_time
from sunpy.net import hek, Fido, attrs as a

###############################################################################
# Let's first grab GOES XRS data for a particular time of interest

tr = TimeRange(['2011-06-07 04:00', '2011-06-07 12:00'])
results = Fido.search(a.Time(tr), a.Instrument('XRS'))
results

###############################################################################
# Then download the data and load it into a TimeSeries

files = Fido.fetch(results)
goes = TimeSeries(files)

###############################################################################
# Next lets grab the HEK data for this time from the NOAA Space Weather
# Prediction Center (SWPC)

client = hek.HEKClient()
flares_hek = client.search(hek.attrs.Time(tr.start, tr.end),
                           hek.attrs.FL, hek.attrs.FRM.Name == 'SWPC')

###############################################################################
# Finally lets plot everything together

fig, ax = plt.subplots()
goes.plot()
ax.axvline(parse_time(flares_hek[0].get('event_peaktime')).plot_date)
ax.axvspan(parse_time(flares_hek[0].get('event_starttime')).plot_date,
           parse_time(flares_hek[0].get('event_endtime')).plot_date,
           alpha=0.2, label=flares_hek[0].get('fl_goescls'))
ax.legend(loc=2)
ax.set_yscale('log')

plt.show()
