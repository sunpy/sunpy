# -*- coding: utf-8 -*-
"""
Test Generic TimeSeries

Created on Thu Jun 23 12:29:55 2016

@author: alex_
"""

import os
import pytest
import datetime
import warnings

import numpy as np
import astropy.units as u
from pandas.util.testing import assert_frame_equal
from collections import OrderedDict

import sunpy
import sunpy.timeseries
from sunpy.time import TimeRange

testpath = sunpy.data.test.rootdir


@pytest.fixture
def eve_test_ts():
    #ToDo: return sunpy.timeseries.TimeSeries(os.path.join(testpath, filename), source='EVE')
    return sunpy.timeseries.TimeSeries(sunpy.data.sample.EVE_LIGHTCURVE, source='EVE')

@pytest.fixture
def fermi_gbm_test_ts():
    #ToDo: return sunpy.timeseries.TimeSeries(os.path.join(testpath, filename), source='GBMSummary')
    return sunpy.timeseries.TimeSeries(sunpy.data.sample.GBM_LIGHTCURVE, source='GBMSummary')

@pytest.fixture
def norrh_test_ts():
    #ToDo: return sunpy.timeseries.TimeSeries(os.path.join(testpath, filename), source='NoRH')
    return sunpy.timeseries.TimeSeries(sunpy.data.sample.NORH_LIGHTCURVE, source='NoRH')

@pytest.fixture
def goes_test_ts():
    #ToDo: return sunpy.timeseries.TimeSeries(os.path.join(testpath, filename), source='GOES')
    return sunpy.timeseries.TimeSeries(sunpy.data.sample.GOES_LIGHTCURVE, source='GOES')

@pytest.fixture
def lyra_test_ts():
    #ToDo: return sunpy.timeseries.TimeSeries(os.path.join(testpath, filename), source='LYRA')
    return sunpy.timeseries.TimeSeries(sunpy.data.sample.LYRA_LEVEL3_LIGHTCURVE, source='LYRA')

@pytest.fixture
def rhessi_test_ts():
    #ToDo: return sunpy.timeseries.TimeSeries(os.path.join(testpath, filename), source='RHESSI')
    return sunpy.timeseries.TimeSeries(sunpy.data.sample.RHESSI_LIGHTCURVE, source='RHESSI')

@pytest.fixture
def noaa_ind_test_ts():
    #ToDo: return sunpy.timeseries.TimeSeries(os.path.join(testpath, filename), source='NOAAPredictIndices')
    return sunpy.timeseries.TimeSeries(sunpy.data.sample.NOAAINDICES_LIGHTCURVE, source='NOAAPredictIndices')

@pytest.fixture
def noaa_pre_test_ts():
    #ToDo: return sunpy.timeseries.TimeSeries(os.path.join(testpath, filename), source='NOAAIndices')
    return sunpy.timeseries.TimeSeries(sunpy.data.sample.NOAAPREDICT_LIGHTCURVE, source='NOAAIndices')

@pytest.fixture
def generic_test_ts():
    #ToDo: generate a generic TS!
    #########
    return sunpy.timeseries.TimeSeries(sunpy.data.sample.NOAAPREDICT_LIGHTCURVE, source='NOAAIndices')


#==============================================================================
# Test TimeSeries Parameters
#==============================================================================

def test_units():
    assert isinstance(eve_test_ts.units, OrderedDict)
    assert isinstance(fermi_gbm_test_ts.units, OrderedDict)
    assert isinstance(norrh_test_ts.units, OrderedDict)
    assert isinstance(goes_test_ts.units, OrderedDict)
    assert isinstance(lyra_test_ts.units, OrderedDict)
    assert isinstance(rhessi_test_ts.units, OrderedDict)
    assert isinstance(noaa_ind_test_ts.units, OrderedDict)
    assert isinstance(noaa_pre_test_ts.units, OrderedDict)
    assert isinstance(generic_test_ts.units, OrderedDict)







    # ToDo: check length? (should match the number of columns)
    # ToDo: is this a good way of going?

#==============================================================================
# Test Truncation Operations
#==============================================================================

@pytest.fixture
def truncation_slice_test_ts_1(eve_test_ts):
    # Truncate by slicing the second half off.
    return eve_test_ts.truncate(0, 720, None)

@pytest.fixture
def truncation_slice_test_ts_2(eve_test_ts):
    # Truncate by slicing the first half off.
    return eve_test_ts.truncate(720, None, None)

def test_truncation_slices(eve_test_ts, truncation_slice_test_ts_1, truncation_slice_test_ts_2):
    # Test resulting DataFrame are similar
    assert len(eve_test_ts.data) == (len(truncation_slice_test_ts_1.data) + len(truncation_slice_test_ts_2.data))
    # Test column lists and unit dictionaries match
    assert eve_test_ts.columns == truncation_slice_test_ts_1.columns == truncation_slice_test_ts_2.columns
    assert eve_test_ts.meta.columns == truncation_slice_test_ts_1.meta.columns == truncation_slice_test_ts_2.meta.columns
    assert eve_test_ts.units == truncation_slice_test_ts_1.units == truncation_slice_test_ts_2.units
    # Test MetaDict match
    assert eve_test_ts.meta.metadata[0][2] == truncation_slice_test_ts_1.meta.metadata[0][2] == truncation_slice_test_ts_2.meta.metadata[0][2]
    # For TS and meta, Test time ranges match for the start and end of the TS.
    assert truncation_slice_test_ts_1.time_range.start == truncation_slice_test_ts_1.meta.timerange.start == eve_test_ts.start
    assert truncation_slice_test_ts_2.time_range.end == truncation_slice_test_ts_2.meta.timerange.end == eve_test_ts.end


@pytest.fixture
def truncation_timerange_test_ts(eve_test_ts):
    # Truncate using a TimeRange object.
    return eve_test_ts.truncate(TimeRange('2012-06-20 02:00:00','2012-06-20 4:00:00'))

def test_truncation_timerange(truncation_timerange_test_ts):
    # Check the resulting timerange in both TS and TSMD
    assert truncation_timerange_test_ts.time_range.hours == truncation_timerange_test_ts.meta.timerange.hours == u.Quantity(2, u.h)

@pytest.fixture
def truncation_dates_test_ts(eve_test_ts):
    # Truncate using strings for start and end datetime.
    return eve_test_ts.truncate('2012-06-20 02:00:00','2012-06-20 4:00:00')

def test_truncation_dates(truncation_dates_test_ts):
    # Check the resulting timerange in both TS and TSMD
    assert truncation_dates_test_ts.time_range.hours == truncation_dates_test_ts.meta.timerange.hours == u.Quantity(2, u.h)


#==============================================================================
# Test Concatenation Operations
#==============================================================================

@pytest.fixture
def concatenated_slices_test_ts(truncation_slice_test_ts_1, truncation_slice_test_ts_2):
    # Concatenate the slices to make a TS similar to the original
    return truncation_slice_test_ts_1.concatenate(truncation_slice_test_ts_2)

def test_concatenation_of_slices(eve_test_ts, concatenated_slices_test_ts):
    # Test resulting DataFrame is similar to the original
    assert_frame_equal(concatenated_slices_test_ts.data, eve_test_ts.data)
    # Otherwise: concatenated_ts.data.equals(eve_test_ts)
    # Compare timeranges from before and after match for both metadata and TS
    assert eve_test_ts.meta.timerange == concatenated_slices_test_ts.meta.timerange
    assert eve_test_ts.time_range == concatenated_slices_test_ts.time_range
    # Test metadata MetaDict matches
    eve_test_ts.meta.metadata[0][2] == concatenated_slices_test_ts.meta.metadata[0][2] == concatenated_slices_test_ts.meta.metadata[1][2]
    # ToDo: Will TSMD.concatenate() want to re-merge the metadata entries back into one?


# Or asseperate tests:
#def test_concatenation_of_slices_dataframe(eve_test_ts, concatenated_slices_test_ts):
#def test_concatenation_of_slices_timerange(eve_test_ts, concatenated_slices_test_ts):
#def test_concatenation_of_slices_metadata(eve_test_ts, concatenated_slices_test_ts):

@pytest.fixture
def concatenation_different_data_test_ts(eve_test_ts, fermi_gbm_test_ts):
    # Tate two different data sources and concatenate
    return eve_test_ts.concatenate(fermi_gbm_test_ts)

def test_concatenation_of_different_data(eve_test_ts, fermi_gbm_test_ts, concatenation_different_data_test_ts):
    # Test data
    # ToDo: figure out some simple tests
    # Test metadata concatenation
    assert concatenation_different_data_test_ts.meta[0] == eve_test_ts.meta[0]
    assert concatenation_different_data_test_ts.meta[1] == fermi_gbm_test_ts.meta[0]
    # Test unit concatenation
    # ToDo: check the resulting units is the union of the originals.




    # ToDo: check units dict is the concatenation of both.


#==============================================================================
# Test Resample Operations
#==============================================================================


    # Subsampling: every n'th element
    ts_eve_subsam = ts_eve.truncate(0,50000,10)
    ts_gbm_subsam = ts_gbm.truncate(0,50000,10)
    ts_norrh_subsam = ts_norrh.truncate(0,50000,10)
    ts_goes_subsam = ts_goes.truncate(0,50000,10)
    ts_lyra_subsam = ts_lyra.truncate(0,50000,10)
    ts_rhessi_subsam = ts_rhessi.truncate(0,50000,10)

    # Resampling: taking an existing lightcurve and resample it at different times to get a new lightcurve
    # summing every 'n' elements of the original time-series
    resampled = ts_eve.resample('3T')         # '3T = 3 mins'
    resampled = ts_eve.resample('H', 'sum') # '60S = 60 seconds'
    resampled = ts_eve.resample('H', 'std')   # 'H = one hour'
    resampled = ts_eve.resample('D', 'mean')  # 'D = one day'
    resampled = ts_eve.resample('D', 'other')  # should show a warning
    # Note: use of the methods: .sum(), .mean(), .std()
    # Note: not sure how to do every nth element.
    #

    # Upsampling: increasing the cadence with interpolation.
    upsampled = ts_eve_subsam.resample('30S', 'pad')
    upsampled = ts_eve_subsam.resample('30S', 'bfill')
    upsampled = ts_eve_subsam.resample('30S', 'ffill')








#==============================================================================
# Test Rotation WCS conversion
#==============================================================================