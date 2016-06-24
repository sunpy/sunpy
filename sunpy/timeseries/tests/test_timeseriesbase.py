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

import sunpy
import sunpy.timeseries
from pandas.util.testing import assert_frame_equal

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

    




#==============================================================================
# Test Truncation Operations
#==============================================================================

@pytest.fixture
def truncation_slice_test_ts_1(eve_test_ts):
    # Alternate values, starting from the first. Derived by slicing.
    return eve_test_ts.truncate(0,100000,2)

@pytest.fixture
def truncation_slice_test_ts_2(eve_test_ts):
    # Alternate values, starting from the second. Derived by slicing.
    return eve_test_ts.truncate(1,100000,2)

    # Truncation: getting shorter duration lightcurves from an existing timeseries. (exclusive of end)
    trunc_str_dates = gts.truncate('20150419','20160607 01')

    tr = TimeRange('2015/04/19', '2016/06/07')
    trunc_time_range = gts.truncate(tr)

    import datetime
    st = datetime.datetime(2016, 6, 7, 0, 0)
    en = datetime.datetime(2016, 6, 8, 12, 0)
    trunc_datetimes = gts.truncate(st,en)

    trunc_slice = gts.truncate(10,500)

    # Subsampling: taking specific elements of an existing lightcurve (for example every n'th element), and creating a new lightcurve
    sub1 = gts.truncate(0,5000,3)




#==============================================================================
# Test Concatenation Operations
#==============================================================================

def test_concatenation_slice_ts(eve_test_ts, truncation_slice_test_ts_1, truncation_slice_test_ts_2):
    # Try to recreate the original TimeSeries, using the two sliced variants.
    concatenated_ts = truncation_slice_test_ts_1.concatenate(truncation_slice_test_ts_2)
    
    # Compare the original ts to the one created using concatenation.
    assert_frame_equal(concatenated_ts, eve_test_ts)
    """
    try:
        assert_frame_equal(eve_test_ts.data, copy.data)
        print('True')
    except:  # appeantly AssertionError doesn't catch all
        print('False')
    """
    # ToDo: check units dict is correctly formed?
    
def concatenation_different_data_test_ts(eve_test_ts, fermi_gbm_test_ts):
    # Tate two different data sources and concatenate
    concatenated_ts = eve_test_ts.concatenate(fermi_gbm_test_ts)
    
    # Now compare using some simple tests.
    assert concatenated_ts.data.shape[0] == eve_test_ts.data.shape[0] + fermi_gbm_test_ts.data.shape[0]
    assert concatenated_ts.data.shape[1] == eve_test_ts.data.shape[1] + fermi_gbm_test_ts.data.shape[1]
    
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





#==============================================================================
# Test Rotation WCS conversion
#==============================================================================