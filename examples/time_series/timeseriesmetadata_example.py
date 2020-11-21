"""
============================
The TimeSeriesMetaData class
============================

This provides an overview of the basic functionality of the TimeSeriesMetaData class.
"""
import astropy.units as u

import sunpy.data.sample
import sunpy.timeseries
from sunpy.net import Fido
from sunpy.net import attrs as a
from sunpy.time import TimeRange, parse_time

##############################################################################
# Search for Data

goes_res = Fido.search(a.Time("2010-11-03", "2010-11-04"), a.Instrument.xrs, a.goes.SatelliteNumber(15))
norh_res = Fido.search(a.Time("2010-11-03", "2010-11-04"), a.Instrument.norh,
                       a.Wavelength(17 * u.GHz))
##############################################################################
# Download the data and load it
goes_files = Fido.fetch(goes_res)
norh_files = Fido.fetch(norh_res)

goes_ts = sunpy.timeseries.TimeSeries(goes_files, source='XRS', concatenate=True)
norh_ts = sunpy.timeseries.TimeSeries(norh_files, source='NoRH', concatenate=True)

##############################################################################
# Combining the two series
large_ts = goes_ts.concatenate(norh_ts)

##############################################################################
# The metadata can be easily viewed:
large_ts.meta

##############################################################################
# You can reduce the depth of the view:
print(large_ts.meta.to_string(2))

##############################################################################
# The TimeSeriesMetaData class stores all the individual file metadata MetaDict
# objects as 3-tuple entries in it's internal list with the TimeRange, list of
# column names and metadictionary. This list is stored in order of ascending
# TR.start.
# Access of the the list is done using the metadata property:
large_ts.meta.metadata

##############################################################################
# The TimeSeriesMetaData class has a number of other properties, including the
# timerange property that returns a TimeRange for the entire metadata:

large_ts.meta.time_range
# Further properties can be used to get lists of details, such as:
# lists of the values:
large_ts.meta.timeranges  # List of the time ranges
large_ts.meta.columns  # List of the column names
large_ts.meta.metas  # List of the meta dictionaries

##############################################################################
# When you truncate the TimeSeries, the metadata is truncated too:
large_trunc_ts = large_ts.truncate(TimeRange('2010-11-03 13:59:57.468999',
                                             '2010-11-04 13:59:56.091999'))
print(large_trunc_ts.meta.to_string(2))

##############################################################################
# Finding metadata can be achieved using the find method and applying filters for
# time and/or colname. This returns another TimeSeriesMetaData object:

large_trunc_ts.meta.find(time=parse_time('2010-11-04 09:01:16'))
large_trunc_ts.meta.find(time='2010-11-04 09:01:16', colname='xrsb')

# You can get the time of a row a from the TimeSeries object's index:
large_trunc_ts.meta.find(time=large_trunc_ts.index[10])

# Note: with no filters you get a duplicate of the original TimeSeriesMetaData
# object.

##############################################################################
# There is also a get method:
large_trunc_ts.meta.get('telescop')

# Again, filters can be used:
large_trunc_ts.meta.get('telescop', time='2010-11-04 09:01:16', colname='xrsb')

# And if we just want the values, the values method returns just a list:
large_trunc_ts.meta.get('telescop').values()

##############################################################################
# You can update values similar to dictionaries, though all of the contained
# MetaDict objects will be updated that match your filter criteria:

large_trunc_ts.meta.update({'new_key_1': 'added to all.'})
large_trunc_ts.meta.update({'new_key_2': 'added to some.'}, colname='xrsa')
print(large_trunc_ts.meta.to_string(2))

# but you can't overwrite previous entries without setting the overwrite kwarg,
# this is to protect the integrity of the metadata:

large_trunc_ts.meta.update({'new_key_1': 'changed'}, overwrite=True)
print(large_trunc_ts.meta.to_string(2))
