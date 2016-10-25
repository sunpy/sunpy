# -*- coding: utf-8 -*-
"""
Created on Wed Jul 20 10:24:06 2016

"""


from __future__ import print_function, division

import copy

from sunpy.timeseries import TimeSeriesMetaData
from sunpy.time import TimeRange
from sunpy.util.metadata import MetaDict
from collections import OrderedDict

import pytest


#==============================================================================
# Creating TimeSeriesMetaData Objects
#==============================================================================

@pytest.fixture
def basic_1_md():
    tr = TimeRange('2010-01-01 13:59:57.468999', '2010-01-02 13:59:56.091999')
    colnames = [ 'column1', 'column2' ]
    metadict = MetaDict(OrderedDict([('md1_key1', 'value1'), ('md1_key2', 'value2'), ('all_same', 'value3'), ('all_different', 'diff_1')]))
    lis = [ ( tr, colnames, metadict ) ]
    return TimeSeriesMetaData(lis)

@pytest.fixture
def basic_2_md():
    tr = TimeRange('2010-01-02 13:59:57.468999', '2010-01-03 13:59:56.091999')
    colnames = [ 'column1', 'column2' ]
    metadict = MetaDict(OrderedDict([('md2_key1', 'value1'), ('md2_key2', 'value2'), ('all_same', 'value3'), ('all_different', 'diff_2')]))
    lis = [ ( tr, colnames, metadict ) ]
    return TimeSeriesMetaData(lis)

@pytest.fixture
def basic_3_md():
    tr = TimeRange('2010-01-03 13:59:57.468999', '2010-01-03 13:59:56.091999')
    colnames = [ 'column1', 'column2' ]
    metadict = MetaDict(OrderedDict([('md3_key1', 'value1'), ('md3_key2', 'value2'), ('all_same', 'value3'), ('all_different', 'diff_3')]))
    lis = [ ( tr, colnames, metadict ) ]
    return TimeSeriesMetaData(lis)

@pytest.fixture
def basic_4_md():
    tr = TimeRange('2010-01-01 20:59:57.468999', '2010-01-03 20:59:56.091999')
    colnames = [ 'md4_column1', 'md4_column2' ]
    metadict = MetaDict(OrderedDict([('md4_key1', 'value1'), ('md4_key2', 'value2'), ('all_same', 'value3'), ('all_different', 'diff_4')]))
    tup = ( tr, colnames, metadict )
    return TimeSeriesMetaData(tup)


@pytest.fixture
def overlap_and_interleave_with_basic_1_md():
    tr = TimeRange('2010-01-01 01:01:00.0', '2010-01-02 01:01:00.0')
    colnames = [ 'column1', 'column2' ]
    metadict = MetaDict(OrderedDict([('other_key1', 'value1'), ('other_key2', 'value2'), ('all_same', 'value3'), ('all_different', 'diff_5')]))
    lis = [ ( tr, colnames, metadict ) ]
    return TimeSeriesMetaData(lis)

#==============================================================================
# Test Creating TimeSeriesMetaData With Limited Input
#==============================================================================

def test_create_mithout_metadata():
    tr = TimeRange('2010-01-01 13:59:57.468999', '2010-01-02 13:59:56.091999')
    colnames = [ 'column1', 'column2' ]
    tsmd_1 = TimeSeriesMetaData(timerange=tr, colnames=colnames)
    assert isinstance(tsmd_1, TimeSeriesMetaData)
    assert tsmd_1.metadata[0][1] == colnames
    tsmd_2 = TimeSeriesMetaData(timerange=tr)
    assert isinstance(tsmd_1, TimeSeriesMetaData)
    assert tsmd_2.metadata[0][1] == []
    assert tsmd_1.metadata[0][0] == tsmd_2.metadata[0][0] == tr
    assert tsmd_1.metadata[0][2] == tsmd_2.metadata[0][2] == MetaDict()
    assert len(tsmd_1.metadata) == len(tsmd_2.metadata) == 1

def test_create_mithout_metadata_or_timerange():
    # without a timerange we should get errors
    colnames = [ 'column1', 'column2' ]
    with pytest.raises(ValueError):
        TimeSeriesMetaData(colnames=colnames)
    with pytest.raises(ValueError):
        TimeSeriesMetaData()


#==============================================================================
# Test Appending TimeSeriesMetaData Objects
#==============================================================================

def test_append_similar(basic_1_md):
    appended = copy.deepcopy(basic_1_md)
    appended.append(*basic_1_md.metadata[0])
    # The duplicate should not have been added
    assert appended == basic_1_md

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
    assert basic_ascending_append_md.metadata[2] == basic_3_md.metadata[0]

def test_basic_descending_append_md(basic_1_md, basic_2_md, basic_3_md, basic_ascending_append_md):
    appended = copy.deepcopy(basic_3_md)
    appended.append(*basic_1_md.metadata[0])
    appended.append(*basic_2_md.metadata[0])
    assert appended == basic_ascending_append_md


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

def test_append_invalid_timerange(basic_1_md):
    appended = copy.deepcopy(basic_1_md)
    with pytest.raises(ValueError):
        appended.append('not_a_timerange', basic_1_md.metadata[0][1], basic_1_md.metadata[0][2])


#==============================================================================
# Test TimeSeriesMetaData Concatenation
#==============================================================================

def test_concatenate(basic_ascending_append_md, basic_4_md, complex_append_md):
    concatenated = copy.deepcopy(basic_ascending_append_md)
    concatenated = concatenated.concatenate(basic_4_md)
    assert concatenated == complex_append_md


#==============================================================================
# Test TimeSeriesMetaData Truncation
#==============================================================================

@pytest.fixture
def truncated_none_md(basic_ascending_append_md):
    # This timerange covers the whole range of metadata, so no change is expected
    tr = TimeRange('2010-01-01 1:59:57.468999', '2010-01-04 23:59:56.091999')
    truncated = copy.deepcopy(basic_ascending_append_md)
    truncated._truncate(tr)
    return truncated

@pytest.fixture
def truncated_start_md(basic_ascending_append_md):
    # This time range starts after the original, so expect truncation
    tr = TimeRange('2010-01-02 20:59:57.468999', '2010-01-04 23:59:56.091999')
    truncated = copy.deepcopy(basic_ascending_append_md)
    truncated._truncate(tr)
    return truncated

@pytest.fixture
def truncated_end_md(basic_ascending_append_md):
    # This time range ends before the original, so expect truncation
    tr = TimeRange('2010-01-01 1:59:57.468999', '2010-01-03 1:59:56.091999')
    truncated = copy.deepcopy(basic_ascending_append_md)
    truncated._truncate(tr)
    return truncated

@pytest.fixture
def truncated_both_md(basic_ascending_append_md):
    # This time range starts after and ends before the original, so expect truncation
    tr = TimeRange('2010-01-02 20:59:57.468999', '2010-01-03 1:59:56.091999')
    truncated = copy.deepcopy(basic_ascending_append_md)
    truncated._truncate(tr)
    return truncated

@pytest.fixture
def truncated_new_tr_all_before_md(basic_ascending_append_md):
    # Time range begins and ends before the data
    tr = TimeRange('2010-01-01 01:01:01.000000', '2010-01-01 02:01:01.000000')
    truncated = copy.deepcopy(basic_ascending_append_md)
    truncated._truncate(tr)
    return truncated

@pytest.fixture
def truncated_new_tr_all_after_md(basic_ascending_append_md):
    # Time range begins and ends after the data
    tr = TimeRange('2010-01-04 01:01:01.000000', '2010-01-04 02:01:01.000000')
    truncated = copy.deepcopy(basic_ascending_append_md)
    truncated._truncate(tr)
    return truncated


#==============================================================================
# Test TimeSeriesMetaData TimeRanges
#==============================================================================

def test_truncated_none_tr(basic_ascending_append_md, truncated_none_md):
    assert basic_ascending_append_md.time_range == truncated_none_md.time_range

def test_truncated_start_tr(truncated_start_md):
    tr = TimeRange('2010-01-02 20:59:57.468999', truncated_start_md.time_range.end)
    assert truncated_start_md.time_range == tr

def test_truncated_end_tr(basic_ascending_append_md, truncated_end_md):
    tr = TimeRange(truncated_end_md.time_range.start, '2010-01-03 1:59:56.091999')
    assert truncated_end_md.time_range == tr

def test_truncated_both_tr(truncated_both_md):
    tr = TimeRange('2010-01-02 20:59:57.468999', '2010-01-03 1:59:56.091999')
    assert truncated_both_md.time_range == tr

def test_truncated_tr_outside(truncated_new_tr_all_before_md, truncated_new_tr_all_after_md):
    assert truncated_new_tr_all_before_md.metadata == truncated_new_tr_all_after_md.metadata == []

def test_basic_ascending_append_tr(basic_1_md, basic_3_md, basic_ascending_append_md):
    tr = TimeRange(basic_1_md.time_range.start, basic_3_md.time_range.end)
    assert basic_ascending_append_md.time_range == tr

def test_complex_append_tr(basic_1_md, basic_4_md, complex_append_md):
    tr = TimeRange(basic_1_md.time_range.start, basic_4_md.time_range.end)
    assert complex_append_md.time_range == tr


#==============================================================================
# Test TimeSeriesMetaData find method
#==============================================================================

def test_find_return_type(basic_ascending_append_md):
    assert isinstance(basic_ascending_append_md.find(), TimeSeriesMetaData)

def test_find_no_filters(basic_ascending_append_md, basic_1_md):
    assert basic_ascending_append_md.find() == basic_ascending_append_md

def test_find_time_filter(basic_1_md, basic_2_md, basic_4_md, complex_append_md):
    assert complex_append_md.find(time='2010-01-01 14:59:57.468999') == basic_1_md
    temp_md = copy.deepcopy(basic_2_md)
    temp_md = temp_md.concatenate(basic_4_md)
    assert complex_append_md.find(time='2010-01-02 20:59:57.468999') == temp_md

def test_find_colname_filter(basic_4_md, complex_append_md, basic_ascending_append_md):
    assert complex_append_md.find(colname='md4_column2') == basic_4_md
    assert complex_append_md.find(colname='column1') == basic_ascending_append_md

def test_find_both_filters(basic_2_md, basic_4_md, complex_append_md):
    assert complex_append_md.find(time='2010-01-02 20:59:57.468999', colname='column2') == basic_2_md
    assert complex_append_md.find(time='2010-01-02 20:59:57.468999', colname='md4_column1') == basic_4_md


#==============================================================================
# Test TimeSeriesMetaData get and update methods
#==============================================================================

def test_get_return_type(complex_append_md):
    assert isinstance(complex_append_md.get('md1_key1'),TimeSeriesMetaData)

def test_get_no_filters(complex_append_md):
    assert complex_append_md.get('md1_key1').values() == ['value1']
    assert complex_append_md.get('all_same').values() == ['value3']

def test_get_time_filter(complex_append_md):
    assert complex_append_md.get('md1_key1', time='2010-01-01 20:59:57.468999').values() == ['value1']
    assert complex_append_md.get('md2_key2', time='2010-01-02 20:59:57.468999').values() == ['value2']
    assert complex_append_md.get('all_same', time='2010-01-01 20:59:57.468999').values() == ['value3']
    assert complex_append_md.get('all_different', time='2010-01-01 20:59:57.468999').values() == ['diff_1', 'diff_4']

def test_get_colname_filter(complex_append_md):
    assert complex_append_md.get('md1_key1', colname='column1').values() == ['value1']
    assert complex_append_md.get('md2_key2', colname='column2').values() == ['value2']
    assert complex_append_md.get('all_same', colname='column1').values() == ['value3']
    assert complex_append_md.get('all_different', colname='column1').values() == ['diff_1', 'diff_2', 'diff_3']

def test_get_both_filters(complex_append_md):
    assert complex_append_md.get('all_different', time='2010-01-02 20:59:57.468999', colname='column2').values() == ['diff_2']


#==============================================================================
# Test TimeSeriesMetaData update method
#==============================================================================

def test_update(complex_append_md):
    assert isinstance(complex_append_md.get('md1_key1'),TimeSeriesMetaData)

def test_update_all_dict_type_input(complex_append_md):
    # Check all three dictionary types work the same
    updated_dict_md = copy.deepcopy(complex_append_md)
    updated_dict_md.update({'added':'added'})
    updated_ordereddict_md = copy.deepcopy(complex_append_md)
    updated_ordereddict_md.update(MetaDict(OrderedDict([('added', 'added')])))
    updated_metadict_md = copy.deepcopy(complex_append_md)
    updated_metadict_md.update(MetaDict(OrderedDict([('added', 'added')])))
    assert updated_dict_md == updated_ordereddict_md == updated_metadict_md

def test_update_overwrite(complex_append_md, overwrite=True):
    updated_not_overwritten_md = copy.deepcopy(complex_append_md)
    updated_not_overwritten_md.update({'all_same': 'updated'})
    updated_overwritten_md = copy.deepcopy(complex_append_md)
    updated_overwritten_md.update({'all_same': 'updated'}, overwrite=True)
    assert updated_not_overwritten_md == complex_append_md
    assert updated_overwritten_md.get('all_same').values() == ['updated']

def test_update_time_filter(complex_append_md):
    updated_md = copy.deepcopy(complex_append_md)
    updated_md.update({'all_same': 'updated'}, time='2010-01-01 20:59:57.468999', overwrite=True)
    assert updated_md.metadata[0][2]['all_same'] == updated_md.metadata[1][2]['all_same'] == 'updated'
    assert updated_md.metadata[2][2]['all_same'] == updated_md.metadata[3][2]['all_same'] == 'value3'

def test_update_colname_filter(complex_append_md):
    updated_md = copy.deepcopy(complex_append_md)
    updated_md.update({'all_same': 'updated'}, colname='column1', overwrite=True)
    assert updated_md.metadata[0][2]['all_same'] == updated_md.metadata[2][2]['all_same'] == updated_md.metadata[3][2]['all_same'] == 'updated'
    assert updated_md.metadata[1][2]['all_same'] == 'value3'

def test_update_both_filters(complex_append_md):
    updated_md = copy.deepcopy(complex_append_md)
    updated_md.update({'all_same': 'updated'}, colname='column1', time='2010-01-01 23:59:57.468999', overwrite=True)
    assert updated_md.metadata[0][2]['all_same'] == 'updated'
    assert updated_md.metadata[1][2]['all_same'] == updated_md.metadata[2][2]['all_same'] == updated_md.metadata[3][2]['all_same'] == 'value3'


#==============================================================================
# Test Misc Methods
#==============================================================================

def test_get_index(basic_1_md, basic_ascending_append_md):
    assert basic_ascending_append_md.get_index(0) == basic_1_md.get_index(0)

def test_equality(basic_1_md, basic_2_md, basic_ascending_append_md):
    basic_1_copy_md = copy.deepcopy(basic_1_md)
    assert basic_1_md == basic_1_copy_md
    assert basic_1_md != basic_2_md
    assert basic_1_md != basic_ascending_append_md

def test_to_string_basic(basic_1_md):
    default_str = basic_1_md.to_string()
    assert isinstance(default_str, str)

    # check this matches the __str__ and __repr__ methods
    default_str == basic_1_md.__str__() == basic_1_md.__repr__()

def test_to_string_depth(basic_1_md):
    depth_1_str = basic_1_md.to_string(depth=1)
    assert len(depth_1_str.split('\n')) == 7
    assert len(depth_1_str.split('...')) == 4
    assert basic_1_md.metadata[0][1][0] in depth_1_str

    # for depth of 2 the time-range will take exactly 2 lines
    depth_2_str = basic_1_md.to_string(depth=2)
    assert len(depth_2_str.split('\n')) == 8
    assert len(depth_2_str.split('...')) == 2
    assert (basic_1_md.metadata[0][1][0] in depth_2_str) and (basic_1_md.metadata[0][1][1] in depth_2_str)

def test_to_string_width(basic_1_md):
    width_110_str = basic_1_md.to_string(width=110)
    split = width_110_str.split('\n')
    assert len(split[0]) == len(split[1]) == len(split[2]) == len(split[3]) == 110

    width_60_str = basic_1_md.to_string(width=60)
    split = width_60_str.split('\n')
    assert len(split[0]) == len(split[1]) == len(split[2]) == len(split[3]) == 60

def test_to_string_few_metadict_entries(basic_1_md):
    tr = basic_1_md.metadata[0][0]
    colnames = basic_1_md.metadata[0][1]
    metadict = MetaDict(OrderedDict([('md1_key1', 'value1')]))
    lis = [ ( tr, colnames, metadict ) ]
    basic_1_md_less_metadict_entries = TimeSeriesMetaData(lis)
    depth_3_str = basic_1_md_less_metadict_entries.to_string(depth = 3)
    assert isinstance(depth_3_str, str)

def test_timeranges(basic_ascending_append_md):
    lis = []
    lis.append(basic_ascending_append_md.metadata[0][0])
    lis.append(basic_ascending_append_md.metadata[1][0])
    lis.append(basic_ascending_append_md.metadata[2][0])
    assert basic_ascending_append_md.timeranges == lis

def test_columns(complex_append_md):
    lis = complex_append_md.metadata[0][1] + complex_append_md.metadata[1][1]
    lis.sort()
    assert complex_append_md.columns == lis

def test_metas(complex_append_md):
    lis = []
    lis.append(complex_append_md.metadata[0][2])
    lis.append(complex_append_md.metadata[1][2])
    lis.append(complex_append_md.metadata[2][2])
    lis.append(complex_append_md.metadata[3][2])
    assert complex_append_md.metas == lis

def test_rename_column(complex_append_md):
    col_renamed_md = copy.deepcopy(complex_append_md)
    old = complex_append_md.metadata[0][1][0]
    new = 'renamed'
    col_renamed_md._rename_column(old, new)
    assert col_renamed_md.metadata[0][1][0] == col_renamed_md.metadata[2][1][0] == col_renamed_md.metadata[3][1][0] == new

def test_remove_column(complex_append_md):
    col_removed_md = copy.deepcopy(complex_append_md)
    col = complex_append_md.metadata[0][1][0]
    col_removed_md._remove_columns(col)
    assert col_removed_md.metadata[0][1] == col_removed_md.metadata[2][1] == col_removed_md.metadata[3][1] == [complex_append_md.metadata[0][1][1]]
    assert col_removed_md.metadata[1][1] == [complex_append_md.metadata[1][1][0], complex_append_md.metadata[1][1][1]]

def test_remove_columns(complex_append_md):
    cols_removed_md = copy.deepcopy(complex_append_md)
    lis = [ complex_append_md.metadata[0][1][0], complex_append_md.metadata[1][1][1]]
    cols_removed_md._remove_columns(lis)
    assert cols_removed_md.metadata[0][1] == cols_removed_md.metadata[2][1] == cols_removed_md.metadata[3][1] == [complex_append_md.metadata[0][1][1]]
    assert cols_removed_md.metadata[1][1] == [complex_append_md.metadata[1][1][0]]

def test_validate_meta_good(complex_append_md):
    assert complex_append_md._validate_meta(complex_append_md)

def test_validate_meta_interleaved(basic_1_md, overlap_and_interleave_with_basic_1_md):
    concatenated = copy.deepcopy(basic_1_md)
    concatenated = concatenated.concatenate(overlap_and_interleave_with_basic_1_md)
    assert concatenated._validate_meta(concatenated)