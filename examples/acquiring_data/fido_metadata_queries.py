"""
====================================
Querying Metadata clients using Fido
====================================

This example shows how to search and retrieve metadata using `~sunpy.net.Fido` from the
search facilities like `~sunpy.net.hek.HEKClient`, `~sunpy.net.helio.HECClient`,
and `~sunpy.net.jsoc.JSOCClient`. It also shows how to display desired columns in the result.
"""
from sunpy.net import Fido
from sunpy.net import attrs as a

###################################################################
# We will query Helio for the 'rhessi_flare_list' table and
# shall limit the total number of records to 5.
# For the same time range, we will query HEK for 'FL' as the
# Event Type and 'PeakFlux' greater than 1000.
# We will also search JSOC for 'hmi.m_45s' Series.
timerange = a.Time('2010/8/1 03:40', '2010/8/1 3:40:10')
results = Fido.search(timerange, a.helio.TableName('rhessi_hxr_flare') |
                      a.hek.FL & (a.hek.FL.PeakFlux > 1000) |
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
hecresults = results.get_response('hec')
hekresults = results.get_response('hek')
jsocresults = results.get_response('jsoc')

###################################################################
# "hekresults" has a lot of columns, we can use ``show()``
# to specify the column names to be displayed.
hektable = hekresults.show('event_peaktime', 'obs_instrument', 'fl_peakflux')
print(hektable)

###################################################################
# ``['T_REC', 'TELESCOP', 'INSTRUME', 'WAVELNTH', 'CAR_ROT']`` are
# default columns shown in `~sunpy.net.jsoc.JSOCResponse`.
# To display all columns from JSOC results, we can use ``show()``
# without passings any arguments to the method.
print(jsocresults)
jsoctable = jsocresults.show()
print(jsoctable)
