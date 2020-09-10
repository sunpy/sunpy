"""
====================================
Querying Metadata clients using Fido
====================================

This example shows how to search and retrieve metadata using `~sunpy.net.Fido` from
search facilities like `~sunpy.net.hek.HEKClient`, `~sunpy.net.helio.HECClient`,
and `~sunpy.net.jsoc.JSOCClient`. It also shows how to display desired columns from the result.
"""
from sunpy.net import Fido
from sunpy.net import attrs as a

###################################################################
# We will query Helio for the 'rhessi_flare_list' table and
# limit the total number of records to 5.
# For the same time range, we will query HEK for 'FL' as the
# event type and 'PeakFlux' greater than 1000.
# We will also search JSOC for 'hmi.m_45s' series.
timerange = a.Time('2010/8/1 03:40', '2010/8/1 3:40:10')
results = Fido.search(timerange, a.helio.TableName('rhessi_hxr_flare') |
                      a.hek.FL & (a.hek.FL.PeakFlux > 1000) |
                      a.jsoc.Series('hmi.m_45s') & a.jsoc.Notify("jsoc@cadair.com"))

###################################################################
# ``results`` is a `~sunpy.net.fido_factory.UnifiedResponse` object that
# contains records returned from querying various clients by "Fido.search".
print(results)

###################################################################
# Now we will download the searched records. Since HEK and HELIO
# clients don't provide files, `Fido.fetch` will
# ignore them and only download files from JSOC.
files = Fido.fetch(results)
print(files)

###################################################################
# Now we will extract individual responses from Fido results.
# We can index these results using the client's name (which is case-insensitive).
hec_results, hek_results, jsoc_results = results['hec'], results['hek'], results['jsoc']

###################################################################
# The results from a metadata search could have up to 100 columns.
# As a result, you can use use ``show()`` to specify the column names you want to display.
hek_table = hek_results.show('event_peaktime', 'obs_instrument', 'fl_peakflux')
print(hek_table)

###################################################################
# The results from JSOC have a default set of columns to show and are
# ``['T_REC', 'TELESCOP', 'INSTRUME', 'WAVELNTH', 'CAR_ROT']``.
# To display all of the columns, we can use ``show()``without passings any arguments.
print(jsoc_results)
jsoc_table = jsoc_results.show()
print(jsoc_table)
