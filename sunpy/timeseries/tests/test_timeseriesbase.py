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
import copy

import numpy as np
import astropy.units as u
from pandas.util.testing import assert_frame_equal
from pandas import DataFrame
import pandas as pd
from collections import OrderedDict
from astropy.tests.helper import assert_quantity_allclose
from astropy.table import Table
from astropy.time import Time

import sunpy
import sunpy.timeseries
from sunpy.time import TimeRange
from sunpy.util.metadata import MetaDict
from sunpy.timeseries import TimeSeriesMetaData
from sunpy.tests.helpers import figure_test, skip_wcsaxes

import sunpy.data.sample
import sunpy.data.test
testpath = sunpy.data.test.rootdir


#==============================================================================
# TimeSeries Tests
#==============================================================================

filepath = sunpy.data.test.rootdir

eve_filepath = os.path.join(filepath, 'EVE_L0CS_DIODES_1m_truncated.txt')
fermi_gbm_filepath = os.path.join(filepath, 'gbm.fits')
norrh_filepath = os.path.join(filepath, '')
goes_filepath = os.path.join(filepath, 'goes.fits')
lyra_filepath = os.path.join(filepath, '')
rhessi_filepath = os.path.join(filepath, '')
noaa_ind_filepath = os.path.join(filepath, '')
noaa_pre_filepath = os.path.join(filepath, '')

goes_filepath = sunpy.data.sample.GOES_LIGHTCURVE
norrh_filepath = sunpy.data.sample.NORH_LIGHTCURVE
lyra_filepath = sunpy.data.sample.LYRA_LEVEL3_LIGHTCURVE
rhessi_filepath = sunpy.data.sample.RHESSI_LIGHTCURVE
noaa_ind_filepath = sunpy.data.sample.NOAAINDICES_LIGHTCURVE
noaa_pre_filepath = sunpy.data.sample.NOAAPREDICT_LIGHTCURVE


@pytest.fixture
def eve_test_ts():
    #ToDo: return sunpy.timeseries.TimeSeries(os.path.join(testpath, filename), source='EVE')
    return sunpy.timeseries.TimeSeries(eve_filepath, source='EVE')

@pytest.fixture
def fermi_gbm_test_ts():
    #ToDo: return sunpy.timeseries.TimeSeries(os.path.join(testpath, filename), source='GBMSummary')
    return sunpy.timeseries.TimeSeries(fermi_gbm_filepath, source='GBMSummary')

@pytest.fixture
def norrh_test_ts():
    #ToDo: return sunpy.timeseries.TimeSeries(os.path.join(testpath, filename), source='NoRH')
    return sunpy.timeseries.TimeSeries(norrh_filepath, source='NoRH')

@pytest.fixture
def goes_test_ts():
    #ToDo: return sunpy.timeseries.TimeSeries(os.path.join(testpath, filename), source='GOES')
    return sunpy.timeseries.TimeSeries(goes_filepath, source='GOES')

@pytest.fixture
def lyra_test_ts():
    #ToDo: return sunpy.timeseries.TimeSeries(os.path.join(testpath, filename), source='LYRA')
    return sunpy.timeseries.TimeSeries(lyra_filepath, source='LYRA')

@pytest.fixture
def rhessi_test_ts():
    #ToDo: return sunpy.timeseries.TimeSeries(os.path.join(testpath, filename), source='RHESSI')
    return sunpy.timeseries.TimeSeries(rhessi_filepath, source='RHESSI')

@pytest.fixture
def noaa_ind_test_ts():
    #ToDo: return sunpy.timeseries.TimeSeries(os.path.join(testpath, filename), source='NOAAIndices')
    return sunpy.timeseries.TimeSeries(noaa_ind_filepath, source='NOAAIndices')

@pytest.fixture
def noaa_pre_test_ts():
    #ToDo: return sunpy.timeseries.TimeSeries(os.path.join(testpath, filename), source='NOAAPredictIndices')
    return sunpy.timeseries.TimeSeries(noaa_pre_filepath, source='NOAAPredictIndices')

@pytest.fixture
def generic_ts():
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

@pytest.fixture
def table_ts():
    # Generate the data and the corresponding dates
    base = datetime.datetime.today()
    times = Time([base - datetime.timedelta(minutes=x) for x in range(0, 24 * 60)])
    intensity = u.Quantity(np.sin(np.arange(0, 12 * np.pi, ((12 * np.pi) / (24*60)))), u.W/u.m**2)

    # Create the units and meta objects
    units = OrderedDict([('intensity', u.W/u.m**2)])
    meta = MetaDict({'key':'value'})
    tbl_meta = MetaDict({'t_key':'t_value'})

    # Create a suitable mixin qtable
    table = Table([times, intensity], names=['time', 'intensity'], meta=tbl_meta)
    table.add_index('time')

    # Create TS from dataframe and check
    return sunpy.timeseries.TimeSeries(table, meta, units)

#==============================================================================
# Test Resulting TimeSeries Parameters
#==============================================================================

def test_units_type(eve_test_ts, fermi_gbm_test_ts, norrh_test_ts, goes_test_ts, lyra_test_ts, rhessi_test_ts, noaa_ind_test_ts, noaa_pre_test_ts, generic_ts, table_ts):
    assert isinstance(eve_test_ts.units, OrderedDict)
    assert isinstance(fermi_gbm_test_ts.units, OrderedDict)
    assert isinstance(norrh_test_ts.units, OrderedDict)
    assert isinstance(goes_test_ts.units, OrderedDict)
    assert isinstance(lyra_test_ts.units, OrderedDict)
    assert isinstance(rhessi_test_ts.units, OrderedDict)
    assert isinstance(noaa_ind_test_ts.units, OrderedDict)
    assert isinstance(noaa_pre_test_ts.units, OrderedDict)
    assert isinstance(generic_ts.units, OrderedDict)
    assert isinstance(table_ts.units, OrderedDict)

def test_meta_type(eve_test_ts, fermi_gbm_test_ts, norrh_test_ts, goes_test_ts, lyra_test_ts, rhessi_test_ts, noaa_ind_test_ts, noaa_pre_test_ts, generic_ts, table_ts):
    assert isinstance(eve_test_ts.meta, TimeSeriesMetaData)
    assert isinstance(fermi_gbm_test_ts.meta, TimeSeriesMetaData)
    assert isinstance(norrh_test_ts.meta, TimeSeriesMetaData)
    assert isinstance(goes_test_ts.meta, TimeSeriesMetaData)
    assert isinstance(lyra_test_ts.meta, TimeSeriesMetaData)
    assert isinstance(rhessi_test_ts.meta, TimeSeriesMetaData)
    assert isinstance(noaa_ind_test_ts.meta, TimeSeriesMetaData)
    assert isinstance(noaa_pre_test_ts.meta, TimeSeriesMetaData)
    assert isinstance(generic_ts.meta, TimeSeriesMetaData)
    assert isinstance(table_ts.meta, TimeSeriesMetaData)

def test_data_type(eve_test_ts, fermi_gbm_test_ts, norrh_test_ts, goes_test_ts, lyra_test_ts, rhessi_test_ts, noaa_ind_test_ts, noaa_pre_test_ts, generic_ts, table_ts):
    assert isinstance(eve_test_ts.data, DataFrame)
    assert isinstance(fermi_gbm_test_ts.data, DataFrame)
    assert isinstance(norrh_test_ts.data, DataFrame)
    assert isinstance(goes_test_ts.data, DataFrame)
    assert isinstance(lyra_test_ts.data, DataFrame)
    assert isinstance(rhessi_test_ts.data, DataFrame)
    assert isinstance(noaa_ind_test_ts.data, DataFrame)
    assert isinstance(noaa_pre_test_ts.data, DataFrame)
    assert isinstance(generic_ts.data, DataFrame)
    assert isinstance(table_ts.data, DataFrame)

    # ToDo: check length? (should match the number of columns)

#==============================================================================
# Test Truncation Operations
#==============================================================================

@pytest.fixture
def truncation_slice_test_ts_1(eve_test_ts):
    # Truncate by slicing the second half off.
    return eve_test_ts.truncate(0, int(len(eve_test_ts.data)/2), None)

@pytest.fixture
def truncation_slice_test_ts_2(eve_test_ts):
    # Truncate by slicing the first half off.
    return eve_test_ts.truncate(int(len(eve_test_ts.data)/2), None, None)

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
    assert truncation_slice_test_ts_1.time_range.start == truncation_slice_test_ts_1.meta.time_range.start == eve_test_ts.time_range.start
    assert truncation_slice_test_ts_2.time_range.end == truncation_slice_test_ts_2.meta.time_range.end == eve_test_ts.time_range.end


@pytest.fixture
def truncation_timerange_test_ts(eve_test_ts):
    # Truncate using a TimeRange object.
    return eve_test_ts.truncate(eve_test_ts.time_range.split(3)[1])

def test_truncation_timerange(eve_test_ts, truncation_timerange_test_ts):
    # Check the resulting timerange in both TS and TSMD
    assert truncation_timerange_test_ts.time_range == truncation_timerange_test_ts.meta.time_range == eve_test_ts.time_range.split(3)[1]

@pytest.fixture
def truncation_dates_test_ts(eve_test_ts):
    # Truncate using strings for start and end datetime.
    start_str = str(eve_test_ts.time_range.split(3)[1].start)
    end_str = str(eve_test_ts.time_range.split(3)[1].end)
    return eve_test_ts.truncate(start_str, end_str)

def test_truncation_dates(eve_test_ts, truncation_dates_test_ts):
    # Check the resulting timerange in both TS and TSMD
    assert truncation_dates_test_ts.time_range == truncation_dates_test_ts.meta.time_range == eve_test_ts.time_range.split(3)[1]


#==============================================================================
# Test Extraction Operations
#==============================================================================

@pytest.fixture
def extraction_test_ts(eve_test_ts):
    # Extract the CMLon column
    return eve_test_ts.extract('CMLon')

def test_extraction(eve_test_ts, extraction_test_ts):
    # Test there's only one column in the data, metadata and units
    assert len(extraction_test_ts.data.columns) == 1
    assert len(extraction_test_ts.meta.columns) == 1
    assert len(extraction_test_ts.units) == 1

    # Test this column name matches
    assert eve_test_ts.data.columns[0] == eve_test_ts.data.columns[0] == list(eve_test_ts.units.keys())[0]

    # Test the data matches
    extracted_df = DataFrame(eve_test_ts.data['CMLon']).dropna()
    extracted_df = extracted_df.sort_index()
    assert_frame_equal(extraction_test_ts.data, extracted_df)

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
    assert eve_test_ts.meta.time_range == concatenated_slices_test_ts.meta.time_range
    assert eve_test_ts.time_range == concatenated_slices_test_ts.time_range
    # Test metadata MetaDict matches
    eve_test_ts.meta.metadata[0][2] == concatenated_slices_test_ts.meta.metadata[0][2] == concatenated_slices_test_ts.meta.metadata[1][2]
    # ToDo: Will TSMD.concatenate() want to re-merge the metadata entries back into one?

@pytest.fixture
def concatenation_different_data_test_ts(eve_test_ts, fermi_gbm_test_ts):
    # Take two different data sources and concatenate
    return eve_test_ts.concatenate(fermi_gbm_test_ts)

def test_concatenation_of_different_data(eve_test_ts, fermi_gbm_test_ts, concatenation_different_data_test_ts):
    # ToDo: test the metadata is as expected using the below. (note ATM this fails if the order is changed)
    #assert concatenation_different_data_test_ts.meta.metadata[0] == fermi_gbm_test_ts.meta.metadata[0]
    #assert concatenation_different_data_test_ts.meta.metadata[1] == eve_test_ts.meta.metadata[0]
    value = True
    for key in list(concatenation_different_data_test_ts.meta.metadata[0][2].keys()):
        if concatenation_different_data_test_ts.meta.metadata[0][2][key] != fermi_gbm_test_ts.meta.metadata[0][2][key]:
            value = False
    for key in list(concatenation_different_data_test_ts.meta.metadata[1][2].keys()):
        if concatenation_different_data_test_ts.meta.metadata[1][2][key] != eve_test_ts.meta.metadata[0][2][key]:
            value = False
    assert value

    # Test units concatenation
    comined_units = copy.deepcopy(eve_test_ts.units)
    comined_units.update(fermi_gbm_test_ts.units)
    assert dict(concatenation_different_data_test_ts.units) == dict(comined_units)

    # Test data is the concatenation
    comined_df = pd.concat([eve_test_ts.data, fermi_gbm_test_ts.data])
    comined_df = comined_df.sort_index()
    assert_frame_equal(concatenation_different_data_test_ts.data, comined_df)

def test_concatenation_different_data_error(eve_test_ts, fermi_gbm_test_ts):
    # Take two different data sources and concatenate but set with the same_source
    # kwarg as true, this should not concatenate.
    with pytest.raises(TypeError):
        eve_test_ts.concatenate(fermi_gbm_test_ts, same_source=True)

def test_generic_construction_concatenation():
    # Generate the data and the corrisponding dates
    base = datetime.datetime.today()
    times = [base - datetime.timedelta(minutes=x) for x in range(0, 24 * 60)]
    intensity1 = np.sin(np.arange(0, 12 * np.pi, ((12 * np.pi) / (24*60))))
    intensity2 = np.sin(np.arange(0, 12 * np.pi, ((12 * np.pi) / (24*60))))

    # Create the data DataFrame, header MetaDict and units OrderedDict
    data = DataFrame(intensity1, index=times, columns=['intensity'])
    data2 = DataFrame(intensity2, index=times, columns=['intensity2'])
    units = OrderedDict([('intensity', u.W/u.m**2)])
    units2 = OrderedDict([('intensity', u.W/u.m**2)])
    meta = MetaDict({'key':'value'})
    meta2 = MetaDict({'key2':'value2'})

    # Create TS individually
    ts_1 = sunpy.timeseries.TimeSeries(data, meta, units)
    ts_2 = sunpy.timeseries.TimeSeries(data2, meta2, units2)
    ts_concat = ts_1.concatenate(ts_2, axis=1)
    assert isinstance(ts_concat, sunpy.timeseries.timeseriesbase.GenericTimeSeries)
    assert len(ts_concat.data) == len(times)
    assert ts_concat.columns == ['intensity', 'intensity2']
    assert len(ts_concat.meta.metadata) == 2

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
    assert ((eve_test_ts.data['CMLon'].values == column_quantity.value) | (np.isnan(eve_test_ts.data['CMLon'].values) & np.isnan(column_quantity.value))).all()

@pytest.fixture
def add_column_from_quantity_ts(eve_test_ts, column_quantity):
    # Add a column to a TS using an astropy quantity
    return eve_test_ts.add_column('quantity_added', column_quantity)

def test_add_column_from_quantity(eve_test_ts, add_column_from_quantity_ts, column_quantity):
    # Test the column similar to the original quantity?
    assert_quantity_allclose(add_column_from_quantity_ts.quantity('quantity_added'), column_quantity)
    # Test the full list of columns are pressent
    assert set(add_column_from_quantity_ts.data.columns) == set(eve_test_ts.data.columns) | set(['quantity_added'])

@pytest.fixture
def add_column_from_array_ts(eve_test_ts, column_quantity):
    # Add a column to a TS using a numpy array
    return eve_test_ts.add_column('array_added', column_quantity.value, unit=column_quantity.unit)

def test_add_column_from_array(eve_test_ts, add_column_from_array_ts, column_quantity):
    # Test the column similar to the original quantity?
    assert_quantity_allclose(add_column_from_array_ts.quantity('array_added'), column_quantity)

    # Test the full list of columns are pressent
    assert set(add_column_from_array_ts.data.columns) == set(eve_test_ts.data.columns) | set(['array_added'])

def test_add_column_from_array_no_units(eve_test_ts, column_quantity):
    ts = eve_test_ts.add_column('array_added', column_quantity.value)
    assert (ts.quantity('array_added') == column_quantity.value).all()

#==============================================================================
# Test Exporting to different formats
#==============================================================================

def test_ts_to_table(generic_ts):
    tbl = generic_ts.to_table()
    assert isinstance(tbl, Table)
    assert tbl.keys() == [ 'date', generic_ts.columns[0] ]
    assert len(tbl) == len(generic_ts.data)
    assert (tbl[generic_ts.columns[0]].quantity == generic_ts.quantity(generic_ts.columns[0])).all()

def test_ts_to_dataframe(generic_ts):
    df = generic_ts.to_dataframe()
    assert isinstance(df, DataFrame)
    assert_frame_equal(df, generic_ts.data)

def test_ts_to_array(generic_ts):
    arr = generic_ts.to_array()
    assert isinstance(arr, np.ndarray)
    assert len(arr) == len(generic_ts.data)

#==============================================================================
# Test Other Functions
#==============================================================================

def test_equality(generic_ts, table_ts):
    generic_copy_ts = copy.deepcopy(generic_ts)
    assert generic_ts == generic_copy_ts
    generic_copy_ts.meta.metadata[0][2]['key'] = 1
    assert generic_ts != generic_copy_ts
    assert generic_ts != table_ts

def test_equality_different_ts_types(generic_ts, table_ts):
    # this should fail as they're not the smae type and can't match
    assert not (generic_ts == eve_test_ts)

def test_ts_index(generic_ts):
    assert (generic_ts.index == generic_ts.data.index).all()

def test_ts_sort_index(generic_ts):
    assert generic_ts.sort_index().data.equals(generic_ts.data.sort_index())

#_validate_units

#_validate_meta

@skip_wcsaxes
@figure_test
def test_eve_peek(eve_test_ts):
    eve_test_ts.peek()


# ToDo:
###Extracting column as quantity or array#ts_eve = ts_eve.add_column(colname, qua_new, overwrite=True)
###Updating a column using quantity or array#ts_eve = ts_eve.add_column(colname, qua_new, overwrite=True)
###Updating the units# ts_eve = ts_eve.add_column(colname, qua_new, unit=unit, overwrite=True)
