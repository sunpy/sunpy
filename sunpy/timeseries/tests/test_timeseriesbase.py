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
from pandas import DataFrame
from collections import OrderedDict

import sunpy
import sunpy.timeseries
from sunpy.time import TimeRange
from sunpy.util.metadata import MetaDict
from sunpy.timeseries import TimeSeriesMetaData

import sunpy.data.sample
import sunpy.data.test
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
    #ToDo: return sunpy.timeseries.TimeSeries(os.path.join(testpath, filename), source='NOAAIndices')
    return sunpy.timeseries.TimeSeries(sunpy.data.sample.NOAAINDICES_LIGHTCURVE, source='NOAAIndices')

@pytest.fixture
def noaa_pre_test_ts():
    #ToDo: return sunpy.timeseries.TimeSeries(os.path.join(testpath, filename), source='NOAAPredictIndices')
    return sunpy.timeseries.TimeSeries(sunpy.data.sample.NOAAPREDICT_LIGHTCURVE, source='NOAAPredictIndices')

@pytest.fixture
def generic_test_ts():
    # Generate the data and the corrisponding dates
    base = datetime.datetime.today()
    dates = [base - datetime.timedelta(minutes=x) for x in range(0, 24 * 60)]
    intensity = np.sin(np.arange(0, 12 * np.pi, ((12 * np.pi) / (24*60))))

    # Create the data DataFrame, header MetaDict and units OrderedDict
    data = DataFrame(intensity, index=dates, columns=['intensity'])
    units = OrderedDict([('intensity', u.W/u.m**2)])
    meta = MetaDict({'key':'value'})

    # Create the time series
    return sunpy.timeseries.TimeSeries(data, meta, units)

#==============================================================================
# Test Creating TimeSeries From Various Dataformats
#==============================================================================

# ToDo:
###ts_goes.to_dataframe()
###ts_goes.to_table()
###ts_goes.to_array()

#==============================================================================
# Test TimeSeries Parameters
#==============================================================================

def test_units_type(eve_test_ts, fermi_gbm_test_ts, norrh_test_ts, goes_test_ts, lyra_test_ts, rhessi_test_ts, noaa_ind_test_ts, noaa_pre_test_ts, generic_test_ts):
    assert isinstance(eve_test_ts.units, OrderedDict)
    assert isinstance(fermi_gbm_test_ts.units, OrderedDict)
    assert isinstance(norrh_test_ts.units, OrderedDict)
    assert isinstance(goes_test_ts.units, OrderedDict)
    assert isinstance(lyra_test_ts.units, OrderedDict)
    assert isinstance(rhessi_test_ts.units, OrderedDict)
    assert isinstance(noaa_ind_test_ts.units, OrderedDict)
    assert isinstance(noaa_pre_test_ts.units, OrderedDict)
    assert isinstance(generic_test_ts.units, OrderedDict)

def test_meta_type(eve_test_ts, fermi_gbm_test_ts, norrh_test_ts, goes_test_ts, lyra_test_ts, rhessi_test_ts, noaa_ind_test_ts, noaa_pre_test_ts, generic_test_ts):
    assert isinstance(eve_test_ts.meta, TimeSeriesMetaData)
    assert isinstance(fermi_gbm_test_ts.meta, TimeSeriesMetaData)
    assert isinstance(norrh_test_ts.meta, TimeSeriesMetaData)
    assert isinstance(goes_test_ts.meta, TimeSeriesMetaData)
    assert isinstance(lyra_test_ts.meta, TimeSeriesMetaData)
    assert isinstance(rhessi_test_ts.meta, TimeSeriesMetaData)
    assert isinstance(noaa_ind_test_ts.meta, TimeSeriesMetaData)
    assert isinstance(noaa_pre_test_ts.meta, TimeSeriesMetaData)
    assert isinstance(generic_test_ts.meta, TimeSeriesMetaData)

def test_data_type(eve_test_ts, fermi_gbm_test_ts, norrh_test_ts, goes_test_ts, lyra_test_ts, rhessi_test_ts, noaa_ind_test_ts, noaa_pre_test_ts, generic_test_ts):
    assert isinstance(eve_test_ts.data, DataFrame)
    assert isinstance(fermi_gbm_test_ts.data, DataFrame)
    assert isinstance(norrh_test_ts.data, DataFrame)
    assert isinstance(goes_test_ts.data, DataFrame)
    assert isinstance(lyra_test_ts.data, DataFrame)
    assert isinstance(rhessi_test_ts.data, DataFrame)
    assert isinstance(noaa_ind_test_ts.data, DataFrame)
    assert isinstance(noaa_pre_test_ts.data, DataFrame)
    assert isinstance(generic_test_ts.data, DataFrame)

    # ToDo: check length? (should match the number of columns)
# ToDo: is this a good way of going?
# ToDo: test for all:
###ts_lyra.data
###ts_lyra.meta
###ts_lyra.time_range # Returns a SunPy TimeRange object.
###ts_lyra.name
###ts_lyra.nickname
###ts_lyra.detector
###ts_lyra.instrument
###ts_lyra.observatory

"""
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
# Test Extraction Operations
#==============================================================================

@pytest.fixture
def extraction_test_ts(eve_test_ts):
    # Extract the CMLon column
    return eve_test_ts.extract('CMLon')

def test_extraction(eve_test_ts, extraction_test_ts):
    # Test there's only one column in the data, metadata and units
    assert len(eve_test_ts.data.columns) == 1
    assert len(eve_test_ts.meta.columns) == 1
    assert len(eve_test_ts.units) == 1

    # Test this column name matches
    assert eve_test_ts.data.columns[0] == eve_test_ts.data.columns[0] == list(eve_test_ts.units.keys())[0]

    # Test the data matches
    assert_frame_equal(eve_test_ts.data, DataFrame(eve_test_ts.data['CMLon']))

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
# Test Data Manipulation
#==============================================================================

@pytest.fixture
def column_quantity(eve_test_ts):
    # Get the astropy Quantity of values from a column
    return eve_test_ts.quantity('CMLon')

def test_column_quantity(eve_test_ts, column_quantity):
    # Test the units and values match
    assert eve_test_ts.units['CMLon'] == column_quantity.unit
    assert eve_test_ts.data['CMLon'].values == column_quantity.value

@pytest.fixture
def add_column_from_quantity_ts(eve_test_ts, column_quantity):
    # Add a column to a TS using an astropy quantity
    return eve_test_ts.quantity('quantity_added', column_quantity)

def test_add_column_from_quantity(eve_test_ts, add_column_from_quantity_ts, column_quantity):
    # Test the column similar to the original quantity?
    assert add_column_from_quantity_ts.quantity('quantity_added') == column_quantity
    # Test the full list of columns are pressent
    assert set(add_column_from_quantity_ts.data.columns) == set(eve_test_ts.data.columns) | set('quantity_added')

@pytest.fixture
def add_column_from_array(eve_test_ts, column_quantity):
    # Add a column to a TS using a numpy array
    return eve_test_ts.quantity('array_added', column_quantity.value, unit=column_quantity.unit)

def test_add_column_from_array(eve_test_ts, add_column_from_array_ts, column_quantity):
    # Test the column similar to the original quantity?
    assert add_column_from_array_ts.quantity('array_added') == column_quantity
    # Test the full list of columns are pressent
    assert set(add_column_from_array_ts.data.columns) == set(eve_test_ts.data.columns) | set('array_added')



# ToDo:
###Extracting column as quantity or array#ts_eve = ts_eve.add_column(colname, qua_new, overwrite=True)
###Updating a column using quantity or array#ts_eve = ts_eve.add_column(colname, qua_new, overwrite=True)
###Updating the units# ts_eve = ts_eve.add_column(colname, qua_new, unit=unit, overwrite=True)
"""