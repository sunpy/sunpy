# -*- coding: utf-8 -*-
"""
============================================
Interacting with Data Using SunPy TimeSeries
============================================

This is an early run-through of the basic functionality of the SunPy TimeSeries
class.
This is intended primarily to demonstrate the current interface for discussion
of the final implementation. Much of the code will be changes as the class is
developed.
"""

from __future__ import print_function, division

import os
import copy

import sunpy.data.sample
import sunpy.timeseries
from sunpy.time import TimeRange, parse_time
from sunpy import lightcurve as lc

import astropy.units as u
from astropy.time import Time
from astropy.table import Table

from collections import OrderedDict
import numpy as np
import datetime
from pandas import DataFrame
from sunpy.util.metadata import MetaDict

##############################################################################
# Creating a TimeSeries from a file can be done using the factory.
ts_eve = sunpy.timeseries.TimeSeries(sunpy.data.sample.EVE_LIGHTCURVE, source='EVE')
ts_goes = sunpy.timeseries.TimeSeries(sunpy.data.sample.GOES_LIGHTCURVE, source='XRS')
ts_lyra = sunpy.timeseries.TimeSeries(sunpy.data.sample.LYRA_LEVEL3_LIGHTCURVE, source='LYRA')
ts_noaa_ind = sunpy.timeseries.TimeSeries(sunpy.data.sample.NOAAINDICES_LIGHTCURVE, source='NOAAIndices')
ts_noaa_pre = sunpy.timeseries.TimeSeries(sunpy.data.sample.NOAAPREDICT_LIGHTCURVE, source='NOAAPredictIndices')
ts_norh = sunpy.timeseries.TimeSeries(sunpy.data.sample.NORH_LIGHTCURVE, source='NoRH')
ts_rhessi = sunpy.timeseries.TimeSeries(sunpy.data.sample.RHESSI_LIGHTCURVE, source='RHESSI')
ts_gbm = sunpy.timeseries.TimeSeries(sunpy.data.sample.GBM_LIGHTCURVE, source='GBMSummary')
# Note: for some FITS files a source can be determined implicitly, however it
# is good practice to delcare it explicitly when possible.

##############################################################################
# You can create a list of TimeSeries objects by using multiple files.
# TimeSeries sources down have the facility to download data for a given time
# (this is meant to be un Unidown), so for now we have to use the old LightCurve
# class to download files that aren’t in the sample data.
goes_lc_1 = lc.GOESLightCurve.create(TimeRange('2012/06/01', '2012/06/02'))
goes_lc_2 = lc.GOESLightCurve.create(TimeRange('2012/06/02', '2012/06/03'))
goes_lc_3 = lc.GOESLightCurve.create(TimeRange('2012/06/03', '2012/06/04'))
filepath_1 = os.path.join(sunpy.config.get('downloads', 'download_dir'), 'go1520120601.fits')
filepath_2 = os.path.join(sunpy.config.get('downloads', 'download_dir'), 'go1520120602.fits')
filepath_3 = os.path.join(sunpy.config.get('downloads', 'download_dir'), 'go1520120603.fits')
# Using these new files you get a list:
lis_goes_ts = sunpy.timeseries.TimeSeries(filepath_1, filepath_2, source='XRS')
lis_goes_ts = sunpy.timeseries.TimeSeries(filepath_1, filepath_2, filepath_3, source='XRS')
# Using concatenate=True kwarg you can merge the files into one TimeSeries:
combined_goes_ts = sunpy.timeseries.TimeSeries(filepath_1, filepath_2, source='XRS', concatenate=True)
combined_goes_ts = sunpy.timeseries.TimeSeries(filepath_1, filepath_2, filepath_3, source='XRS', concatenate=True)
combined_goes_ts.peek()
# Note: ATM we only accept TimeSeries of a single class being created together
# with the factory. The issue is that several source filetypes don't contain
# metadata that enables us to reliably implicitly gather the source and ATM the
# source is given as a single keyword argument for simplicity. But you can merge
# different TimeSeries classes using concatenate.
# Debate: are we OK for one source at a time?

##############################################################################
# You can concatenate manually:
combined_goes_ts = lis_goes_ts[0].concatenate(lis_goes_ts[1])
fig = combined_goes_ts.peek()
# Note: peek returns a matplotlib figure object, which can be saved to a file using:
fig.savefig('figure.png')

##############################################################################
# The TimeSeries object has 3 primary data storage components:
# data (pandas.DataFrame): stores the data.
# meta (OrderedDict): stores the metadata (like the Map)
# units (OrderedDict): stores the units for each column, with keys that match
# the name of each column.
# These can be accessed like on the map:
ts_lyra.data
ts_lyra.meta
ts_lyra.units
# There are a couple of other useful properties you can quickly get:
ts_lyra.time_range
ts_lyra.index
ts_lyra.columns
# Further data is available from within the metadata, you can filter out for a
# key using the TimeSeriesMetaData.get() method:
combined_goes_ts.meta.get('telescop')
# Note: this returns a TimeSeriesMetaData object, to get a list of just the
# values for this key use the values property of the metadata:
combined_goes_ts.meta.get('telescop').values()
# Note: this always returns a list because there may be one or more results.

##############################################################################
# The ID used in the data Pandas DataFrame object will be a datetime, as can
# be seen using ts_lyra.index.
# You can access a specific value within the TimeSeries data DataFrame using
# all the normal Pandas methods.
# For example, the row with the index of 2015-01-01 00:02:00.008000:
ts_lyra.data.loc[parse_time('2015-01-01 00:02:00.008000')]
# Pandas will actually parse a string to a datetime automatically if it can:
ts_lyra.data.loc['2015-01-01 00:02:00.008000']
# Pandas includes methods to find the indexes of the max/min values in a dataframe:
lyra_ch1_max_index = ts_lyra.data['CHANNEL1'].idxmax()
lyra_ch1_min_index = ts_lyra.data['CHANNEL1'].idxmin()

##############################################################################
# The TimeSeriesMetaData can be summarised:
combined_goes_ts.meta
print(combined_goes_ts.meta)
print(combined_goes_ts.meta.to_string(2))

##############################################################################
# The TimeSeries objects can be visualised using peek():
ts_goes.peek()
# And you can use subplots:
ts_eve.peek(subplots=True)

##############################################################################
# An individual column can be extracted from a TimeSeries:
ts_eve_extract = ts_eve.extract('CMLon')
# Note: no matter the source type of the original TimeSeries, the extracted
# TimeSeries is always generic.

##############################################################################
# You can truncate a TimeSeries using the truncate() method.
# This can use string datetime arguments, a SunPy TimeRange or integer value
# arguments (similar to slicing, but using function notation).
# Using integers we can get every other entry using:
ts_goes_trunc = ts_goes.truncate(0,100000,2)
# Or using a TimeRange:
tr = TimeRange('2012-06-01 05:00','2012-06-01 06:30')
ts_goes_trunc = ts_goes.truncate(tr)
# Or using strings:
ts_goes_trunc = ts_goes.truncate('2012-06-01 05:00','2012-06-01 06:30')
ts_goes_trunc.peek()
# Note: the strings are parsed using SunPy's string parser.
# Debate: how should we deal with metadata when truncating.

##############################################################################
# You can use Pandas resample method, for example to downsample:
df_downsampled = ts_goes_trunc.data.resample('10T', 'mean')
# To get this into a similar TimeSeries we can copy the original:
ts_downsampled = copy.deepcopy(ts_goes_trunc)
ts_downsampled.data = df_downsampled
ts_downsampled.peek()
# You can use 'mean', 'sum' and 'std' methods and any other methods in Pandas.
# Note: changing values within the datframe directly will often affect the units
# involved, but these won't be picked up by the TimeSeries object. Take care
# when doing this to ensure dimensional consistancy.

##############################################################################
# Similarly, to upsample:
df_upsampled = ts_downsampled.data.resample('1T', 'ffill')
# And this can be made into a TimeSeries using:
ts_upsampled = copy.deepcopy(ts_downsampled)
ts_upsampled.data = df_upsampled
ts_upsampled.peek()
# Note: 'ffill', 'bfill' and 'pad' methods work, and as before others should also.

##############################################################################
# The data from the TimeSeries can be retrieved in a number of formats:
ts_goes.to_dataframe()
ts_goes.to_table()
ts_goes.to_array()
# Note: the array doesn't include the datetime index column.

##############################################################################
# Creating a TimeSeries from scratch can be done in a lot of ways, much like a
# Map.
# Input data can be in the form of a Pandas DataFrame (preferred), an astropy
# Table or a Numpy Array.
# To generate some data and the corresponding dates
base = datetime.datetime.today()
dates = [base - datetime.timedelta(minutes=x) for x in range(0, 24 * 60)]
intensity = np.sin(np.arange(0, 12 * np.pi, ((12 * np.pi) / (24*60))))
# Create the data DataFrame, header MetaDict and units OrderedDict
data = DataFrame(intensity, index=dates, columns=['intensity'])
units = OrderedDict([('intensity', u.W/u.m**2)])
meta = MetaDict({'key':'value'})
# Create the time series
ts_custom = sunpy.timeseries.TimeSeries(data, meta, units)

# A more manual dataset would be a numpy array, which we can creat using:
tm = Time(['2000:002', '2001:345', '2002:345'])
a = [1, 4, 5]
b = [2.0, 5.0, 8.2]
c = ['x', 'y', 'z']
arr = np.stack([tm, a, b, c], axis=1)
# Note: this array needs to have the times in the first column, this can be in
# any form that can be converted using astropy.time.Time().

# We can use the array directly:
ts_from_arr   = sunpy.timeseries.TimeSeries(arr,{})

# We can use this to create a table and even include units:
t = Table([tm, a, b, c], names=('time', 'a', 'b', 'c'), meta={'name': 'table'})
t['b'].unit = 's' # Adding units
ts_from_table = sunpy.timeseries.TimeSeries(t,{})

# If you wanted to make a dataframe from this array then you could use:
df = DataFrame(data=arr[:,1:])
df.index = tm
ts_from_df    = sunpy.timeseries.TimeSeries(df,{})

##############################################################################
# You can optionally add units data, a dictionary matching column heading keys
# to an astropy unit.
units = OrderedDict([('a', u.Unit("ct")),
             ('b', u.Unit("ct")),
             ('c', u.Unit("ct"))])
ts_from_table = sunpy.timeseries.TimeSeries(t,{}, units)
ts_from_df = sunpy.timeseries.TimeSeries(df,{}, units)

##############################################################################
# Changing the units for a column simply requires changing the value:
ts_from_table.units['a'] = u.m

##############################################################################
# Quantities can be extracted from a column using the quantity(col_name) method:
colname = 'CMLat'
qua = ts_eve.quantity(colname)
print(qua)

##############################################################################
# You can add or overwrite a column using the add_column method.
# This method ascepts an astropy quantity and will convert to the intended units
# if necessary.
qua_new = u.Quantity(qua.value * 0.01, ts_eve.units[colname])
print(qua_new)
ts_eve = ts_eve.add_column(colname, qua_new, overwrite=True)
# Otherwise you can also use a numpy array and it assume you're using the original
# units:
arr_new = u.Quantity(qua.value * 0.1, ts_eve.units[colname]).value
ts_eve = ts_eve.add_column(colname, qua_new, overwrite=True)
# Finally, if you want to change the units used, you can specify a new unit for
# the column using the unit keyword:
qua_new = u.Quantity(qua.value * 0.00001, ts_eve.units[colname])
unit = u.W/(u.km**2)
ts_eve = ts_eve.add_column(colname, qua_new, unit=unit, overwrite=True)
