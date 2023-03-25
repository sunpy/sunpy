"""
==================================
Querying the GOES flare event list
==================================

How to retrieve the GOES flare event list, sort and save the results.
"""
from sunpy.net import Fido
from sunpy.net import attrs as a

###################################################################
# We use Fido to to query the HEK catalogue. We define our event type
# as a flare ("FL"). We also set up the start and end times over which
# we will search for flare events. We want the list of events
# that were detected by the GOES X-ray Sensor (XRS) instrument between
# ``tstart`` and ``tend`` and GOES class > M1.0.

event_type = "FL"
tstart = "2013/10/28"
tend = "2013/10/29"
result = Fido.search(a.Time(tstart, tend),
                     a.hek.EventType(event_type),
                     a.hek.FL.GOESCls > "M1.0",
                     a.hek.OBS.Observatory == "GOES")

###################################################################
# The result is returned as a `~sunpy.net.fido_factory.UnifiedResponse`,
# from which we can see a table from one provider is found and returned.

# Here we only show two columns due there being over 100 columns returned normally.
print(result.show("hpc_bbox", "refs"))

# It"s also possible to access the HEK results from the
# `~sunpy.net.fido_factory.UnifiedResponse` by name.
hek_results = result["hek"]

###################################################################
# We can also print the key names that correspond to the HEK
# parameters returned by the query.

# We only print every 10th key to avoid the output being too long.
print(hek_results.colnames[::10])

###################################################################
# The results returned contain a lot of information and we may
# only want to keep some main results such as start time,
# end time, peak time, GOES-class, and active region number.
# This can be done as so:

filtered_results = hek_results["event_starttime", "event_peaktime",
                               "event_endtime", "fl_goescls", "ar_noaanum"]

###################################################################
# ``hek_result`` is already sorted by date, if one wants to sort
# by magnitude, one can do the following:

# Sorting is done by using the flare class from "fl_goescls"
# By converting the flare class to a number using ord()
# and adding the flare strength, we can sort by value
by_magnitude = sorted(filtered_results, key=lambda x: ord(x['fl_goescls'][0]) + float(x['fl_goescls'][1:]), reverse=True)

for flare in by_magnitude:
    print(f"Class {flare['fl_goescls']} occurred on {flare['event_starttime']}")

###################################################################
# These results can then be saved to a CSV file, or any other file
# format that `~astropy.table.Table` supports.

filtered_results.write("october_M1_flares.csv", format="csv")
