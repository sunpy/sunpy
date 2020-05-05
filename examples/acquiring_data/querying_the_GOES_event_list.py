"""
==================================================
Querying the GOES flare event list through the HEK
==================================================

How to retrieve the GOES flare event list through use of
SunPy's Heliophysics Event Knowledgebase (HEK) client.
"""
from sunpy.net import hek

###################################################################
# We first set up the HEK client and define our event type that we
# will search for, which is flare ('FL'). We also set up the start
# and end times over which we will search for flare events.
client = hek.HEKClient()
event_type = 'FL'

tstart = '2013/10/28'
tend = '2013/10/29'

###################################################################
# We then use the client to query the HEK for a list of events
# that were detected by the GOES X-ray Sensor (XRS) instrument between
# `tstart` and `tend`.
result = client.search(hek.attrs.Time(tstart, tend),
                       hek.attrs.EventType(event_type),
                       hek.attrs.OBS.Observatory == 'GOES')

###################################################################
# The result is returned as a `~sunpy.net.hek.hek.HEKTable`.
# We can print the number of flares and inspect the result information.
# We can also print the key-values that correspond to the HEK parameters
# returned in result.
print(len(result))
print(result)
print(result.keys())

###################################################################
# We can also specify what GOES class flares we want to query.
# For example, let's search only for large solar flares with
# a GOES class > M1.0. This can be achieved by using the FL.GOESCls
# attribute of the HEK client:

result_m1 = client.search(hek.attrs.Time(tstart, tend),
                          hek.attrs.EventType(event_type),
                          hek.attrs.FL.GOESCls > 'M1.0',
                          hek.attrs.OBS.Observatory == 'GOES')

print(result_m1)

###################################################################
# The results returned to the `~sunpy.net.hek.hek.HEKTable`
# contain a lot of information and we may only want to keep some main
# results such as start time, end time, peak time, GOES-class, and
# active region number. This can be done as so:
new_table = result_m1[['event_starttime', 'event_peaktime',
                       'event_endtime', 'fl_goescls', 'ar_noaanum', ]]

###################################################################
# These results can then be saved to a CSV file, or any other file
# format that `~astropy.table.Table` supports
new_table.write('october_M1_flares.csv', format='csv')
