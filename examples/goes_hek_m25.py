"""
===============================
GOES Flare and HEK Plot Example
===============================

An example showing how to combine GOES and HEK data
"""

import matplotlib.pyplot as plt
from sunpy.lightcurve import GOESLightCurve
from sunpy.time import TimeRange, parse_time
from sunpy.net import hek

###############################################################################
# Let's first grab GOES XRS data for a particular time of interest

tr = TimeRange(['2011-06-07 04:00', '2011-06-07 12:00'])
goes = GOESLightCurve.create(tr)

###############################################################################
# Next lets grab the HEK data for this time from the NOAA Space Weather Prediction Center (SWPC)

client = hek.HEKClient()
flares_hek = client.search(hek.attrs.Time(tr.start,tr.end), hek.attrs.FL, hek.attrs.FRM.Name == 'SWPC')

###############################################################################
# Finally lets plot everything together

goes.peek()
plt.axvline(parse_time(flares_hek[0].get('event_peaktime')))
plt.axvspan(parse_time(flares_hek[0].get('event_starttime')), parse_time(flares_hek[0].get('event_endtime')),
            alpha=0.2, label=flares_hek[0].get('fl_goescls'))
plt.legend(loc=2)
