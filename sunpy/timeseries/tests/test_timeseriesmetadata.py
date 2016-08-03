# -*- coding: utf-8 -*-
"""
Created on Wed Jul 20 10:24:06 2016

@author: alex_
"""


from __future__ import print_function, division

import os
import copy

import sunpy.data.sample
import sunpy.timeseries
from sunpy.timeseries import TimeSeriesMetaData
from sunpy.time import TimeRange, parse_time
from sunpy.util.metadata import MetaDict
from collections import OrderedDict

from datetime import timedelta

import pytest

#==============================================================================
# Creating TimeSeriesMetaData Objects
#==============================================================================

@pytest.fixture
def basic_1_md():
    tr = TimeRange('2010-01-01 13:59:57.468999', '2010-01-01 13:59:56.091999')
    colnames = [ 'md1_column1', 'md1_column2' ]
    metadict = MetaDict({'md1_key1':'value1', 'md1_key2':'value2', 'all_same':'value3', 'all_different':'diff_1'})
    metadict = MetaDict(OrderedDict([('md1_key1', 'value1'), ('md1_key2', 'value2'), ('all_same', 'value3'), ('all_different', 'diff_1')]))
    lis = [ ( tr, colnames, metadict ) ]
    return TimeSeriesMetaData(lis)

@pytest.fixture
def basic_2_md():
    tr = TimeRange('2010-01-02 13:59:57.468999', '2010-01-02 13:59:56.091999')
    colnames = [ 'md2_column1', 'md2_column2' ]
    metadict = MetaDict({'md2_key1':'value1', 'md2_key2':'value2', 'all_same':'value3', 'all_different':'diff_2'})
    metadict = MetaDict(OrderedDict([('md2_key1', 'value1'), ('md2_key2', 'value2'), ('all_same', 'value3'), ('all_different', 'diff_2')]))
    lis = [ ( tr, colnames, metadict ) ]
    return TimeSeriesMetaData(lis)

@pytest.fixture
def basic_3_md():
    tr = TimeRange('2010-01-03 13:59:57.468999', '2010-01-03 13:59:56.091999')
    colnames = [ 'md3_column1', 'md3_column2' ]
    metadict = MetaDict({'md3_key1':'value1', 'md3_key2':'value2', 'all_same':'value3', 'all_different':'diff_3'})
    metadict = MetaDict(OrderedDict([('md3_key1', 'value1'), ('md3_key2', 'value2'), ('all_same', 'value3'), ('all_different', 'diff_3')]))
    lis = [ ( tr, colnames, metadict ) ]
    return TimeSeriesMetaData(lis)

@pytest.fixture
def basic_4_md():
    tr = TimeRange('2010-01-02 20:59:57.468999', '2010-01-04 13:59:56.091999')
    colnames = [ 'md4_column1', 'md4_column2' ]
    metadict = MetaDict({'md4_key1':'value1', 'md4_key2':'value2', 'all_same':'value3', 'all_different':'diff_4'})
    metadict = MetaDict(OrderedDict([('md4_key1', 'value1'), ('md4_key2', 'value2'), ('all_same', 'value3'), ('all_different', 'diff_4')]))
    lis = [ ( tr, colnames, metadict ) ]
    return TimeSeriesMetaData(lis)

"""
tr = TimeRange('2010-01-01 13:59:57.468999', '2010-01-01 13:59:56.091999')
colnames = [ 'md1_column1', 'md1_column2' ]
metadict = MetaDict({'md1_key1':'value1', 'md1_key2':'value2', 'all_same':'value3', 'all_different':'diff_1'})
metadict = MetaDict(OrderedDict([('md1_key1', 'value1'), ('md1_key2', 'value2'), ('all_same', 'value3'), ('all_different', 'diff_1')]))
lis = [ ( tr, colnames, metadict ) ]
basic_1_md = TimeSeriesMetaData(lis)

tr = TimeRange('2010-01-02 13:59:57.468999', '2010-01-02 13:59:56.091999')
colnames = [ 'md2_column1', 'md2_column2' ]
metadict = MetaDict({'md2_key1':'value1', 'md2_key2':'value2', 'all_same':'value3', 'all_different':'diff_2'})
lis = [ ( tr, colnames, metadict ) ]
basic_2_md = TimeSeriesMetaData(lis)

tr = TimeRange('2010-01-03 13:59:57.468999', '2010-01-03 13:59:56.091999')
colnames = [ 'md3_column1', 'md3_column2' ]
metadict = MetaDict({'md3_key1':'value1', 'md3_key2':'value2', 'all_same':'value3', 'all_different':'diff_3'})
lis = [ ( tr, colnames, metadict ) ]
basic_3_md = TimeSeriesMetaData(lis)

tr = TimeRange('2010-01-02 20:59:57.468999', '2010-01-04 13:59:56.091999')
colnames = [ 'md4_column1', 'md4_column2' ]
metadict = MetaDict({'md4_key1':'value1', 'md4_key2':'value2', 'all_same':'value3', 'all_different':'diff_4'})
lis = [ ( tr, colnames, metadict ) ]
basic_4_md = TimeSeriesMetaData(lis)

appended = copy.deepcopy(basic_1_md)
appended.append(*basic_2_md.metadata[0])
appended.append(*basic_3_md.metadata[0])
basic_ascending_append_md = appended

appended = copy.deepcopy(basic_1_md)
appended.append(*basic_2_md.metadata[0])
appended.append(*basic_3_md.metadata[0])
appended.append(*basic_4_md.metadata[0])
complex_append_md = appended
"""
#==============================================================================
# Test Appending TimeSeriesMetaData Objects
#==============================================================================

@pytest.fixture
def basic_ascending_append_md(basic_1_md, basic_2_md, basic_3_md):
    appended = copy.deepcopy(basic_1_md)
    appended.append(*basic_2_md.metadata[0])
    appended.append(*basic_3_md.metadata[0])
    return appended

def test_basic_ascending_append_md(basic_1_md, basic_2_md, basic_3_md, basic_ascending_append_md):
    # Check all the entries are in the correct order
    assert basic_ascending_append_md.metadata[0] == basic_1_md.metadata[0]
    assert basic_ascending_append_md.metadata[1] == basic_2_md.metadata[0]

    print('\n\n')
    print(basic_ascending_append_md.metadata[2])
    print('\n\n')
    print(basic_3_md.metadata[0])
    print('\n\n')

    assert basic_ascending_append_md.metadata[2] == basic_3_md.metadata[0]

def test_basic_descending_append_md(basic_1_md, basic_2_md, basic_3_md, basic_ascending_append_md):
    appended = copy.deepcopy(basic_3_md)
    appended.append(*basic_1_md.metadata[0])
    appended.append(*basic_2_md.metadata[0])
    assert appended == basic_ascending_append_md

"""
def test_basic_random_append_md(basic_1_md, basic_2_md, basic_3_md, basic_ascending_append_md):
    appended = copy.deepcopy(basic_3_md)
    appended.append(*basic_1_md.metadata[0])
    appended.append(*basic_2_md.metadata[0])
    assert appended == basic_ascending_append_md

@pytest.fixture
def complex_append_md(basic_1_md, basic_2_md, basic_3_md, basic_4_md):
    appended = copy.deepcopy(basic_1_md)
    appended.append(*basic_2_md.metadata[0])
    appended.append(*basic_3_md.metadata[0])
    appended.append(*basic_4_md.metadata[0])
    return appended

def test_complex_append_md(basic_1_md, basic_2_md, basic_3_md, basic_4_md, complex_append_md):
    # Check all the entries are in the correct order
    assert complex_append_md.metadata[0] == basic_1_md.metadata[0]
    assert complex_append_md.metadata[1] == basic_4_md.metadata[0]
    assert complex_append_md.metadata[2] == basic_2_md.metadata[0]
    assert complex_append_md.metadata[3] == basic_3_md.metadata[0]

#==============================================================================
# Test TimeSeriesMetaData Truncation
#==============================================================================

@pytest.fixture
def truncated_none_md(basic_ascending_append_md):
    tr = TimeRange('2010-01-01 1:59:57.468999', '2010-01-03 23:59:56.091999')
    truncated = copy.deepcopy(basic_ascending_append_md)
    return truncated.truncate(tr)

@pytest.fixture
def truncated_start_md(basic_ascending_append_md):
    tr = TimeRange('2010-01-01 20:59:57.468999', '2010-01-03 23:59:56.091999')
    truncated = copy.deepcopy(basic_ascending_append_md)
    return truncated.truncate(tr)

@pytest.fixture
def truncated_end_md(basic_ascending_append_md):
    tr = TimeRange('2010-01-01 1:59:57.468999', '2010-01-03 1:59:56.091999')
    truncated = copy.deepcopy(basic_ascending_append_md)
    return truncated.truncate(tr)

@pytest.fixture
def truncated_both_md(basic_ascending_append_md):
    tr = TimeRange('2010-01-01 20:59:57.468999', '2010-01-03 1:59:56.091999')
    truncated = copy.deepcopy(basic_ascending_append_md)
    return truncated.truncate(tr)


#==============================================================================
# Test TimeSeriesMetaData TimeRanges
#==============================================================================

def test_truncated_none_tr(basic_ascending_append_md, truncated_none_md):
    assert basic_ascending_append_md.time_range == truncated_none_md.time_range

def test_truncated_start_tr(basic_ascending_append_md, truncated_start_md):
    tr = TimeRange('2010-01-01 20:59:57.468999', truncated_start_md.end)
    assert truncated_start_md.time_range == tr

def test_truncated_end_tr(basic_ascending_append_md, truncated_end_md):
    tr = TimeRange(truncated_start_md.start, '2010-01-03 1:59:56.091999')
    assert truncated_end_md.time_range == tr

def test_truncated_both_tr(truncated_both_md):
    tr = TimeRange('2010-01-01 20:59:57.468999', '2010-01-03 1:59:56.091999')
    assert truncated_both_md.time_range == tr

def test_basic_ascending_append_tr(basic_1_md, basic_3_md, basic_ascending_append_md):
    tr = TimeRange(basic_1_md.start, basic_3_md.end)
    assert basic_ascending_append_md.time_range == tr

def test_complex_append_tr(basic_1_md, basic_4_md, complex_append_md):
    tr = TimeRange(basic_1_md.start, basic_4_md.end)
    assert complex_append_md.time_range == tr


#==============================================================================
# Test TimeSeriesMetaData find and update methods
#==============================================================================

def test_find(basic_ascending_append_md):
    assert isinstance(basic_ascending_append_md.get('md1_key1'),TimeSeriesMetaData)

def test_find_no_filters(basic_ascending_append_md, basic_1_md):
    assert basic_ascending_append_md.find('md1_key1') == basic_1_md
    assert basic_ascending_append_md.get('all_same').values == ['value3']

def test_get_time_filter(basic_ascending_append_md):
    assert basic_ascending_append_md.get('md1_key1', time='2010-01-01 20:59:57.468999').values == ['value1']
    assert basic_ascending_append_md.get('md2_key2', time='2010-01-02 20:59:57.468999').values == ['value2']
    assert basic_ascending_append_md.get('all_same', time='2010-01-01 20:59:57.468999').values == ['value3']
    assert basic_ascending_append_md.get('all_different', time='2010-01-01 20:59:57.468999').values == ['diff_1','diff_2', 'diff_3']

def test_get_colname_filter(basic_ascending_append_md):
    assert basic_ascending_append_md.get('md1_key1', colname='md1_column1').values == ['value1']
    assert basic_ascending_append_md.get('md2_key2', colname='md2_column2').values == ['value2']
    assert basic_ascending_append_md.get('all_same', colname='md1_column1').values == ['value3']
    assert basic_ascending_append_md.get('all_different', colname='md1_column1').values == ['diff_1']

def test_get_both_filters(basic_ascending_append_md):
    assert basic_ascending_append_md.get('all_different', time='2010-01-02 20:59:57.468999', colname='md2_column2').values == ['diff_2']

#==============================================================================
# Test TimeSeriesMetaData get and update methods
#==============================================================================

def test_get(basic_ascending_append_md):
    assert isinstance(basic_ascending_append_md.get('md1_key1'),TimeSeriesMetaData)

def test_get_no_filters(basic_ascending_append_md):
    assert basic_ascending_append_md.get('md1_key1').values == ['value1']
    assert basic_ascending_append_md.get('all_same').values == ['value3']

def test_get_time_filter(basic_ascending_append_md):
    assert basic_ascending_append_md.get('md1_key1', time='2010-01-01 20:59:57.468999').values == ['value1']
    assert basic_ascending_append_md.get('md2_key2', time='2010-01-02 20:59:57.468999').values == ['value2']
    assert basic_ascending_append_md.get('all_same', time='2010-01-01 20:59:57.468999').values == ['value3']
    assert basic_ascending_append_md.get('all_different', time='2010-01-01 20:59:57.468999').values == ['diff_1','diff_2', 'diff_3']

def test_get_colname_filter(basic_ascending_append_md):
    assert basic_ascending_append_md.get('md1_key1', colname='md1_column1').values == ['value1']
    assert basic_ascending_append_md.get('md2_key2', colname='md2_column2').values == ['value2']
    assert basic_ascending_append_md.get('all_same', colname='md1_column1').values == ['value3']
    assert basic_ascending_append_md.get('all_different', colname='md1_column1').values == ['diff_1']

def test_get_both_filters(basic_ascending_append_md):
    assert basic_ascending_append_md.get('all_different', time='2010-01-02 20:59:57.468999', colname='md2_column2').values == ['diff_2']
"""
"""
##############################################################################
# Getting GOES data over a satelite change
start_date = "2010-11-02"
start_date = parse_time(start_date)
end_date = start_date
for i in range(0, 6):
    start_date = end_date
    end_date = start_date + timedelta(days=1)
    goes_lc_path = lc.GOESLightCurve.create(TimeRange(start_date.isoformat()[:10], end_date.isoformat()[:10]))
goes_files = [ 'go1420101102.fits', 'go1420101103.fits', 'go1420101104.fits', 'go1520101105.fits', 'go1520101106.fits', 'go1520101107.fits' ]
goes_filepaths = []
for file in goes_files:
    goes_filepaths.append(os.path.join(sunpy.config.get('downloads', 'download_dir'), file))
goes_ts = sunpy.timeseries.TimeSeries(*goes_filepaths, source='GOES', concatenate=True)

# Getting NORH data over same dates
start_date = "2010-11-02"
start_date = parse_time(start_date)
end_date = start_date
for i in range(0, 6):
    start_date = end_date
    end_date = start_date + timedelta(days=1)
    norh_lc_path = lc.NoRHLightCurve.create(start_date.isoformat()[:10])
norh_files = [ 'tca101102', 'tca101103', 'tca101104', 'tca101105', 'tca101106', 'tca101107' ]
norh_filepaths = []
for file in norh_files:
    norh_filepaths.append(os.path.join(sunpy.config.get('downloads', 'download_dir'), file))
norh_ts = sunpy.timeseries.TimeSeries(*norh_filepaths, source='NoRH', concatenate=True)

# Combining
large_ts = goes_ts.concatenate(norh_ts)
# ToDo: Fix: plot doesn't work, it's of type goes TS and so it doesn't plot non-goes data. Should concanate default to GenericTimeSeries???

##############################################################################
# The metadata can be easily viewed:
large_ts.meta

##############################################################################
# You can reduce the depth of the view:
print(large_ts.meta.to_string(2))

##############################################################################
# The TimeSeriesMetaData class stores all the individual file metadata MetaDict
# objects as 3-tuple entries in it's internal list with the TimeRange, list of
# column names and metadictionary. This list is stores in order of ascending
# TR.start.
# Access of the the list is done using the metadata property:
large_ts.meta.metadata

##############################################################################
# The TimeSeriesMetaData class has a number of other properties, including the
# timerange property that returns a TimeRange for the entire metadata:
large_ts.meta.timerange
# Further properties can be used to get lists of details, such as:
# lists of the values:
large_ts.meta.timeranges # List of the time ranges
large_ts.meta.columns    # List of the column names
large_ts.meta.metas      # List of the meta dictionaries

##############################################################################
# When you truncate the TimeSeries, the metadata is truncated too:
large_trunc_ts = large_ts.truncate(TimeRange('2010-11-03 13:59:57.468999', '2010-11-05 13:59:56.091999'))
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
large_trunc_ts.meta.update({'new_key_1':'added to all.'})
large_trunc_ts.meta.update({'new_key_2':'added to some.'}, colname='xrsa')
print(large_trunc_ts.meta.to_string(2))
# but you can't overwrite previous entries without setting the overwrite kwarg,
# this is to protect the integrity of the metadata:
large_trunc_ts.meta.update({'new_key_1':'changed'}, overwrite=True)
print(large_trunc_ts.meta.to_string(2))
"""