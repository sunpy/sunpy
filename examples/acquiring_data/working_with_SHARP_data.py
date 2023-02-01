"""
==========================================
Querying SHARP data and constructing a Map
==========================================

How to query SHARP data through use of
sunpy's jsoc client and construct a map.
"""
import astropy.units as u

import sunpy.map
from sunpy.net import Fido
from sunpy.net import attrs as a

###################################################################
# We use Fido to query the SHARP data. We want data in the time range
# between ``tstart`` and ``tend``, having the time interval for the
# data sampling as ``sampling_rate``. We include the SHARP series
# which we want to query with ``series``, passing the primekey
# using ``primekey_label`` and ``primekey_value``. When JSOC has
# staged our request, a notification will be sent to the ``notify_email``
# address.
# Sometimes, there can be more than one file present for each record.
# ``segment`` is used then to define which files to download.

tstart = "2011-02-10 22:00:00"
tend = "2011-02-10 22:30:00"
sampling_rate = 1*u.hour
series = "hmi.sharp_cea_720s"
primekey_label = "HARPNUM"
primekey_value = 377
notify_email = "sunnycrockett@sunpy.org"
segment = "Bp"

###################################################################
# We plug all these attributes into Fido to conduct the search.
# Once the search is done we can use Fido to fetch the data.

results = Fido.search(a.Time(tstart, tend),
                      a.Sample(sampling_rate),
                      a.jsoc.Series(series),
                      a.jsoc.PrimeKey(primekey_label, primekey_value),
                      a.jsoc.Notify(notify_email),
                      a.jsoc.Segment(segment))

files = Fido.fetch(results)

###################################################################
# Once the data has been queried, we can construct a map with it
# using ``sunpy.map.Map()``

sharp_bp = sunpy.map.Map(files)
