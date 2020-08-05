"""
====================================
Querying Metadata clients using Fido
====================================

How to perform and inspect metadata only queries
using Fido from the HEK and HELIO clients.
"""
from sunpy.net import Fido, attrs as a

###################################################################
# We will query Helio for the 'rhessi_flare_list' table and
# shall limit the total number of records to 5.
# For the same time range, we will query HEK for 'FL' as the
# Event Type and 'PeakFlux' greater than 1000.
# We will also search JSOC for 'hmi.m_45s' Series.
timerange = a.Time('2010/8/1 01:00', '2010/8/1 18:00')
results = Fido.search(timerange, a.helio.TableName('rhessi_hxr_flare') |
                      a.hek.FL & (a.hek.FL.PeakFlux>1000) |
                      a.jsoc.Series('hmi.m_45s') & a.jsoc.Notify("sunpy@sunpy.org"))
print(results)

###################################################################
# Know we will download the searched records. Since HEK and HELIO
# clients don't need to fetch any file records, so fetch will
# ignore them and only download files for JSOC Client.
files = Fido.fetch(results)
print(files)

###################################################################
# Now we will extract individual responses from Fido results.
hecresults = results.get_response(0)
hekresults = results.get_response(1)
jsocresults = results.get_response(2)

###################################################################
# "hekresults" has a lot of columns, we can use ``show()``
# to specify the column names to be displayed.
hektable = hekresults.show('event_peaktime', 'obs_instrument', 'fl_peakflux')
print(hektable)

###################################################################
# If no arguments are specified in show method, then all columns
# are shown.
jsoctable = jsocresults.show()
print(jsoctable)
