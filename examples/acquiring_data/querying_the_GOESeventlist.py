# -*- coding: utf-8 -*-
"""
==================================================
Querying the GOES flare event list through the HEK
==================================================

How to retrieve the GOES falre event list throguh use of 
SunPy's Heliophysics Event Knowledgebase (HEK) client.
"""

from sunpy.net import hek

client = hek.HEKClient()
event_type = 'FL'


tstart = '2013-10-28'
tend = '2013-10-29'

result=client.search(hek.attrs.Time(tstart,tend),
                     hek.attrs.EventType(event_type),
                     hek.attrs.OBS.Observatory == 'GOES')

###################################################################
# The result is returned as an astropy.table.Table. Lets print the 
# number of flares and inspect the result information. We can also print 
# the key-values that correspond to the HEK parameters returned in 
# result.
print(len(result))
print(result)
print(result.keys())


###################################################################
# We can also specify when GOES class flares we want to query. 
# For example, let's search only for large solar flares with
# a GOES class > M1.0. This can be achieved by using the FL.GOESCls
# attribute of HEK:

goes_class_filter = 'M1.0'

result_m1 = client.search(hek.attrs.Time(tstart, tend),
                       hek.attrs.EventType(event_type),
                       hek.attrs.FL.GOESCls > goes_class_filter,
                       hek.attrs.OBS.Observatory == 'GOES')

print(result_m1)

###################################################################
# The results returned to the HEKTable contain a lot of information 
# and we may only want to keep some main results such as start time, 
# end time, peak time, GOES-class, active region number.

new_table = result_m1[['event_starttime', 'event_peaktime', 
                       'event_endtime', 'fl_goescls', 'ar_noaanum', ]]

###################################################################
# These results can then be saved to a csv file, or whatever file
# format that astropy.table.Table supports. 

new_table.save('october_M1_flares.csv', format='csv')





