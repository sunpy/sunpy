"""
==================================================
Querying the GOES flare event list through the HEK
==================================================

How to retrieve the GOES flare event list through use of
SunPy's Heliophysics Event Knowledgebase (HEK) client.
"""
from sunpy.net import Fido
from sunpy.net import attrs as a

###################################################################
# We use Fido to to query the HEK catalogue. We define our event type
# as a flare ('FL'). We also set up the start and end times over which
# we will search for flare events. We want the list of events
# that were detected by the GOES X-ray Sensor (XRS) instrument between
# `tstart` and `tend`.
event_type = 'FL'
tstart = '2013/10/28'
tend = '2013/10/29'
result = Fido.search(a.Time(tstart, tend),
                     a.hek.EventType(event_type),
                     a.hek.OBS.Observatory == 'GOES')

###################################################################
# The result is returned as a `~sunpy.net.fido_factory.UnifiedResponse`,
# from which we can see a table from one provider is found and returned.
print(result)

# We can now inspect the results in `result` by indexing to return get the data
# from the `~sunpy.net.hek.hek.HEKResponse`. Lets print the number of flares and
# inspect the result information.
print(len(result[0]))
# We can also print the key-values that correspond to the HEK parameters returned
# in result[0]. The .table attribute returns an `~astropy.table.Table`.
print(result[0].table.keys())

###################################################################
# We can also specify what GOES class flares we want to query.
# For example, let's search only for large solar flares with
# a GOES class > M1.0. This can be achieved by using the FL.GOESCls
# attribute of the HEK client:

result_m1 = Fido.search(a.Time(tstart, tend),
                        a.hek.EventType(event_type),
                        a.hek.FL.GOESCls > 'M1.0',
                        a.hek.OBS.Observatory == 'GOES')

print(result_m1[0])

###################################################################
# The results returned to the `~sunpy.net.hek.hek.HEKResponse`
# contain a lot of information and we may only want to keep some main
# results such as start time, end time, peak time, GOES-class, and
# active region number. This can be done as so:
new_table = result_m1[0].table['event_starttime', 'event_peaktime',
                               'event_endtime', 'fl_goescls', 'ar_noaanum']

###################################################################
# These results can then be saved to a CSV file, or any other file
# format that `~astropy.table.Table` supports
new_table.write('october_M1_flares.csv', format='csv')
