# -*- coding: utf-8 -*-
"""
=========================================
Interacting with Data Using SunPy TimeSeries
=========================================

This is an early run-through of the basic functionality of the SunPy TimeSeries
class.
This is intended primarily to demonstrate the current interface for discussion
of the final implementation. Much of the code will be changes as the class is
developed.
"""

from __future__ import print_function, division

import os

import sunpy.data.sample
import sunpy.timeseries
from sunpy.time import TimeRange

from sunpy import lightcurve as lc
import astropy.units as u
from collections import OrderedDict
from astropy.table import Table
import numpy as np

##############################################################################
# Creating a TimeSeries from a file can be done using the factory.
ts_eve = sunpy.timeseries.TimeSeries(sunpy.data.sample.EVE_LIGHTCURVE, source='EVE')
ts_gbm = sunpy.timeseries.TimeSeries(sunpy.data.sample.GBM_LIGHTCURVE, source='GBMSummary')
ts_goes = sunpy.timeseries.TimeSeries(sunpy.data.sample.GOES_LIGHTCURVE, source='GOES')
ts_lyra = sunpy.timeseries.TimeSeries(sunpy.data.sample.LYRA_LEVEL3_LIGHTCURVE, source='LYRA')
ts_noaa_ind = sunpy.timeseries.TimeSeries(sunpy.data.sample.NOAAINDICES_LIGHTCURVE, source='NOAAIndices')
ts_noaa_pre = sunpy.timeseries.TimeSeries(sunpy.data.sample.NOAAPREDICT_LIGHTCURVE, source='NOAAPredictIndices')
ts_norh = sunpy.timeseries.TimeSeries(sunpy.data.sample.NORH_LIGHTCURVE, source='NoRH')
ts_rhessi = sunpy.timeseries.TimeSeries(sunpy.data.sample.RHESSI_LIGHTCURVE, source='RHESSI')
# Note: currently you need to define a source, in future this may become implicit for some sources, depending on if the file contains such data.
# Debate: would it be better for consistency to simply always demand the source?

##############################################################################
# You can create a list of TimeSeries objects by using multiple files.
# TimeSeries sources down have the facility to download data for a given time
# (this is meant to be un Unidown), so for now we have to use the old LightCurve
# class to download files that arnt in the sample data.
goes_lc_1 = lc.GOESLightCurve.create(TimeRange('2012/06/01', '2012/06/02'))
goes_lc_2 = lc.GOESLightCurve.create(TimeRange('2012/06/02', '2012/06/03'))
goes_lc_3 = lc.GOESLightCurve.create(TimeRange('2012/06/03', '2012/06/04'))
filepath_1 = os.path.join(sunpy.config.get('downloads', 'download_dir'), 'go1520120601.fits')
filepath_2 = os.path.join(sunpy.config.get('downloads', 'download_dir'), 'go1520120602.fits')
filepath_3 = os.path.join(sunpy.config.get('downloads', 'download_dir'), 'go1520120603.fits')
# Using these new files you get a list:
lis_goes_ts = sunpy.timeseries.TimeSeries(filepath_1, filepath_2, source='GOES')
lis_goes_ts = sunpy.timeseries.TimeSeries(filepath_1, filepath_2, filepath_3, source='GOES')
# Using concatenate=True kwarg you can merge the files into one TimeSeries:
combined_goes_ts = sunpy.timeseries.TimeSeries(filepath_1, filepath_2, source='GOES', concatenate=True)
combined_goes_ts = sunpy.timeseries.TimeSeries(filepath_1, filepath_2, filepath_3, source='GOES', concatenate=True)
combined_goes_ts.peek()
# Note: ATM we only accept TimeSeries of a single class being created together
# with the factory. The issue is that several source filetypes don't contain
# metadata that enables us to reliably implicitly gather the source and ATM the
# source is given as a single keyword argument for simplicity. But you can merge
# different Timeseries classes using concatenate.
# Debate: are we OK for one source at a time?

##############################################################################
# You can concatenate manually:
combined_goes_ts = lis_goes_ts[0].concatenate(lis_goes_ts[1])
combined_goes_ts.peek()
# Debate: how should we deal with metadata when concatenating.

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
# Further data is avalible:
ts_lyra.time_range # Returns a SunPy TimeRange object.
ts_lyra.name
ts_lyra.nickname
ts_lyra.date
ts_lyra.detector
ts_lyra.dsun
ts_lyra.exposure_time
ts_lyra.instrument
ts_lyra.measurement
ts_lyra.observatory
# Note: much of the behaviour of these is yet to be sorted for specific instruments.

##############################################################################
# The TimeSeries objects can be visualised using peek():
ts_goes.peek()
# And you can use subplots:
ts_eve.peek(subplots=True)

##############################################################################
# An individual column can be extracted from a TimeSeries:
ts_eve_extract = ts_eve.extract('CMLon')

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
# You can call Pandas resample method, for example to downsample:
downsampled = ts_goes_trunc.resample('10T', 'mean')
downsampled.peek()
# Note: you can use 'mean', 'sum' and 'std' methods, and it should work with
# other methods that exist in Pandas.

##############################################################################
# To upsample:
upsampled = downsampled.resample('1T', 'ffill')
upsampled.peek()
# Note: 'ffill', 'bfill' and 'pad' methods work, and as before others should also.

##############################################################################
# The data from the TimeSeries can be retrieved in a number of formats:
ts_goes.to_dataframe()
ts_goes.to_table()
ts_goes.to_array()
# Note: the array doesn't include the datetime index column.

##############################################################################
# Creating a TimeSeries from scratch can be done in a lot of ways, much like a
# Map. Input data can be in the form of a Pandas DataFrame (prefered), an astropy
# Table or a Numpy Array.
arr = np.array([[1,2],[3,4]])
a = [1, 4, 5]
b = [2.0, 5.0, 8.2]
c = ['x', 'y', 'z']
t = Table([a, b, c], names=('a', 'b', 'c'), meta={'name': 'first table'})
t['b'].unit = 's' # Adding units
df = t.to_pandas()
ts_from_arr = sunpy.timeseries.TimeSeries(arr,{})
ts_from_table = sunpy.timeseries.TimeSeries(t,{})
ts_from_df = sunpy.timeseries.TimeSeries(df,{})

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