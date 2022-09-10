"""
========================
Getting data from CDAWeb
========================

How to download data from the Coordinated Data Analysis Web (CDAWeb).

CDAWeb stores data from from current and past space physics missions, and is
full of heliospheric insitu datasets.
"""

from sunpy.net import Fido
from sunpy.net import attrs as a
from sunpy.timeseries import TimeSeries

###############################################################################
# `sunpy.net.Fido` is the primary interface to search for and download data and
# will automatically search CDAWeb when the ``cdaweb.Dataset`` attribute is provided to
# the search. To lookup the different dataset IDs available, you can use the
# form at https://cdaweb.gsfc.nasa.gov/index.html/
trange = a.Time('2021/07/01', '2021/07/08')
dataset = a.cdaweb.Dataset('SOLO_L2_MAG-RTN-NORMAL-1-MINUTE')
result = Fido.search(trange, dataset)

###############################################################################
# Let's inspect the results. We can see that there's seven files, one for each
# day within the query.
print(result)

###############################################################################
# Let's download the first two files
downloaded_files = Fido.fetch(result[0, 0:2])
print(downloaded_files)

###############################################################################
# Finally we can load and take a look at the data using
# `~sunpy.timeseries.TimeSeries` This requires an installation of the cdflib
# Python library to read the CDF file.
solo_mag = TimeSeries(downloaded_files, concatenate=True)
print(solo_mag.columns)
solo_mag.peek(columns=['B_RTN_0', 'B_RTN_1', 'B_RTN_2'])
