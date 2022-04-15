import os
import copy
import glob
import datetime
from collections import OrderedDict

import numpy as np
import pandas as pd
import pytest
from erfa.core import ErfaWarning
from pandas import DataFrame
from pandas.testing import assert_frame_equal

import astropy.units as u
from astropy.table import Table
from astropy.tests.helper import assert_quantity_allclose
from astropy.time import TimeDelta

import sunpy
import sunpy.data.test
import sunpy.timeseries
from sunpy.tests.helpers import figure_test
from sunpy.time import TimeRange, parse_time
from sunpy.timeseries import TimeSeriesMetaData
from sunpy.util import SunpyUserWarning
from sunpy.util.metadata import MetaDict

# =============================================================================
# TimeSeries Tests
# =============================================================================

filepath = sunpy.data.test.rootdir

eve_filepath = os.path.join(filepath, 'EVE_L0CS_DIODES_1m_truncated.txt')
esp_filepath = os.path.join(filepath, 'eve_l1_esp_2011046_00_truncated.fits')
fermi_gbm_filepath = os.path.join(filepath, 'gbm.fits')
norh_filepath = os.path.join(filepath, 'tca110810_truncated')
goes_filepath = os.path.join(filepath, 'goes.fits')
lyra_filepath = os.path.join(filepath,
                             'lyra_20150101-000000_lev3_std_truncated.fits.gz')
rhessi_filepath = os.path.join(filepath,
                               'hsi_obssumm_20120601_018_truncated.fits.gz')
noaa_ind_json_filepath = os.path.join(filepath, 'observed-solar-cycle-indices-truncated.json')
noaa_pre_json_filepath = os.path.join(filepath,
                                      'predicted-solar-cycle-truncated.json')
goes_filepath = os.path.join(filepath, 'go1520120601.fits.gz')

a_list_of_many = glob.glob(os.path.join(filepath, "eve", "*"))


@pytest.fixture
def eve_test_ts():
    with pytest.warns(SunpyUserWarning, match='Unknown units'):
        return sunpy.timeseries.TimeSeries(eve_filepath, source='EVE')


@pytest.fixture
def esp_test_ts():
    return sunpy.timeseries.TimeSeries(esp_filepath, source='ESP')


@pytest.fixture
def fermi_gbm_test_ts():
    return sunpy.timeseries.TimeSeries(fermi_gbm_filepath, source='GBMSummary')


@pytest.fixture
def norh_test_ts():
    return sunpy.timeseries.TimeSeries(norh_filepath, source='NoRH')


@pytest.fixture
def goes_test_ts():
    return sunpy.timeseries.TimeSeries(goes_filepath, source='XRS')


@pytest.fixture
def lyra_test_ts():
    return sunpy.timeseries.TimeSeries(lyra_filepath, source='LYRA')


@pytest.fixture
def rhessi_test_ts():
    return sunpy.timeseries.TimeSeries(rhessi_filepath, source='RHESSI')


@pytest.fixture
def noaa_ind_json_test_ts():
    return sunpy.timeseries.TimeSeries(noaa_ind_json_filepath, source='NOAAIndices')


@pytest.fixture
def noaa_pre_json_test_ts():
    # NOAA pre data contains years long into the future, which ERFA complains about
    with pytest.warns(ErfaWarning, match=r'.*dubious year'):
        return sunpy.timeseries.TimeSeries(
            noaa_pre_json_filepath, source='NOAAPredictIndices')


@pytest.fixture
def generic_ts():
    # Generate the data and the corresponding dates
    base = parse_time("2016/10/01T05:00:00")
    dates = base - TimeDelta(np.arange(24 * 60)*u.minute)
    intensity = np.sin(np.arange(0, 12 * np.pi, ((12 * np.pi) / (24 * 60))))
    intensity2 = np.cos(np.arange(0, 12 * np.pi, ((12 * np.pi) / (24 * 60))))

    # Create the data DataFrame, header MetaDict and units dict
    data = DataFrame(np.column_stack([intensity, intensity2]),
                     index=dates.isot.astype('datetime64'),
                     columns=['intensity', 'intensity2'])
    units = {'intensity': u.W / u.m**2,
             'intensity2': u.W / u.m**2}
    meta = MetaDict({'key': 'value'})

    # Create the time series
    return sunpy.timeseries.TimeSeries(data, meta, units)


@pytest.fixture
def concatenate_multi_files_ts():
    with pytest.warns(SunpyUserWarning, match='Unknown units'):
        return sunpy.timeseries.TimeSeries(
            a_list_of_many, source='EVE', concatenate=True)

# =============================================================================
# Test Creating TimeSeries From Various Dataformats
# =============================================================================


@pytest.fixture
def table_ts():
    # Generate the data and the corresponding dates
    base = parse_time(datetime.datetime.today())
    times = base - TimeDelta(np.arange(24 * 60)*u.minute)
    intensity = u.Quantity(
        np.sin(np.arange(0, 12 * np.pi, ((12 * np.pi) / (24 * 60)))), u.W / u.m ** 2)

    # Create the units and meta objects
    units = OrderedDict([('intensity', u.W / u.m**2)])
    meta = MetaDict({'key': 'value'})
    tbl_meta = MetaDict({'t_key': 't_value'})

    # Create a suitable mixin qtable
    table = Table(
        [times, intensity], names=['time', 'intensity'], meta=tbl_meta)
    table.add_index('time')

    # Create TS from dataframe and check
    return sunpy.timeseries.TimeSeries(table, meta, units)


# =============================================================================
# Test Resulting TimeSeries Parameters
# =============================================================================
@pytest.fixture(params=['eve_test_ts', 'esp_test_ts', 'fermi_gbm_test_ts', 'norh_test_ts', 'goes_test_ts',
                        'lyra_test_ts', 'rhessi_test_ts', 'noaa_ind_json_test_ts',
                        'noaa_pre_json_test_ts', 'generic_ts', 'table_ts'])
def many_ts(request):
    # Fixture to return lots of different timeseries
    return request.getfixturevalue(request.param)


def test_units_type(many_ts):
    assert isinstance(many_ts.units, OrderedDict)


def test_meta_type(many_ts):
    assert isinstance(many_ts.meta, TimeSeriesMetaData)


def test_data_type(many_ts):
    assert isinstance(many_ts.to_dataframe(), DataFrame)
    # TODO: check length? (should match the number of columns)

# =============================================================================
# Test Basic Single-Timeseries Truncation Operations
# =============================================================================


@pytest.fixture
def truncation_slice_test_ts_1(eve_test_ts):
    # Truncate by keeping only the first quarter.
    return eve_test_ts.truncate(0, int(len(eve_test_ts.to_dataframe()) / 4), None)


@pytest.fixture
def truncation_slice_test_ts_2(eve_test_ts):
    # Truncate by keeping only the second quarter.
    return eve_test_ts.truncate(
        int(len(eve_test_ts.to_dataframe()) / 4), int(len(eve_test_ts.to_dataframe()) / 2), None
    )


@pytest.fixture
def truncation_slice_test_ts_3(eve_test_ts):
    # Truncate by keeping only the third quarter.
    return eve_test_ts.truncate(
        int(len(eve_test_ts.to_dataframe()) / 2), int(3 * len(eve_test_ts.to_dataframe()) / 4), None
    )


@pytest.fixture
def truncation_slice_test_ts_4(eve_test_ts):
    # Truncate by keeping only the fourth quarter.
    return eve_test_ts.truncate(
        int(3 * len(eve_test_ts.to_dataframe()) / 4), len(eve_test_ts.to_dataframe()), None
    )


def test_truncation_slices(eve_test_ts,
                           truncation_slice_test_ts_1, truncation_slice_test_ts_2,
                           truncation_slice_test_ts_3, truncation_slice_test_ts_4):
    # Test resulting DataFrame are similar
    assert len(eve_test_ts.to_dataframe()) == (len(truncation_slice_test_ts_1.to_dataframe()) +
                                               len(truncation_slice_test_ts_2.to_dataframe()) +
                                               len(truncation_slice_test_ts_3.to_dataframe()) +
                                               len(truncation_slice_test_ts_4.to_dataframe()))
    # Test column lists and unit dictionaries match
    assert (eve_test_ts.columns ==
            truncation_slice_test_ts_1.columns ==
            truncation_slice_test_ts_2.columns ==
            truncation_slice_test_ts_3.columns ==
            truncation_slice_test_ts_4.columns)
    assert (eve_test_ts.meta.columns ==
            truncation_slice_test_ts_1.meta.columns ==
            truncation_slice_test_ts_2.meta.columns ==
            truncation_slice_test_ts_3.meta.columns ==
            truncation_slice_test_ts_4.meta.columns)
    assert (eve_test_ts.units ==
            truncation_slice_test_ts_1.units ==
            truncation_slice_test_ts_2.units ==
            truncation_slice_test_ts_3.units ==
            truncation_slice_test_ts_4.units)
    # Test MetaDict match
    assert (eve_test_ts.meta.metadata[0][2] ==
            truncation_slice_test_ts_1.meta.metadata[0][2] ==
            truncation_slice_test_ts_2.meta.metadata[0][2] ==
            truncation_slice_test_ts_3.meta.metadata[0][2] ==
            truncation_slice_test_ts_4.meta.metadata[0][2])
    # For TS and meta, Test time ranges match for the start and end of the TS.
    assert (truncation_slice_test_ts_1.time_range.start ==
            truncation_slice_test_ts_1.meta.time_range.start ==
            eve_test_ts.time_range.start)
    assert (truncation_slice_test_ts_4.time_range.end ==
            truncation_slice_test_ts_4.meta.time_range.end ==
            eve_test_ts.time_range.end)


def test_truncation_timerange(eve_test_ts):
    truncated = eve_test_ts.truncate(eve_test_ts.time_range.split(3)[1])
    # Check the resulting timerange in both TS and TSMD
    assert (truncated.time_range ==
            truncated.meta.time_range ==
            eve_test_ts.time_range.split(3)[1])


def test_truncation_dates(eve_test_ts):
    start_str = str(eve_test_ts.time_range.split(3)[1].start)
    end_str = str(eve_test_ts.time_range.split(3)[1].end)
    truncated = eve_test_ts.truncate(start_str, end_str)
    # Check the resulting timerange in both TS and TSMD
    assert (truncated.time_range ==
            truncated.meta.time_range ==
            eve_test_ts.time_range.split(3)[1])


# =============================================================================
# Test Basic Single-Timeseries Truncation Operations
# =============================================================================
def test_truncated_none_ts(concatenate_multi_files_ts):
    # This timerange covers the whole range of metadata, so no change is expected
    a = concatenate_multi_files_ts.meta.metadata[0][0].start - TimeDelta(1*u.day)
    b = concatenate_multi_files_ts.meta.metadata[-1][0].end + TimeDelta(1*u.day)
    tr = TimeRange(a, b)
    truncated = copy.deepcopy(concatenate_multi_files_ts)
    truncated = truncated.truncate(tr)

    assert concatenate_multi_files_ts.meta == truncated.meta


def test_truncated_start_ts(concatenate_multi_files_ts):
    a = concatenate_multi_files_ts.meta.metadata[1][0].center
    b = concatenate_multi_files_ts.meta.metadata[-1][0].end + TimeDelta(1*u.day)
    tr = TimeRange(a, b)
    truncated = copy.deepcopy(concatenate_multi_files_ts)
    truncated = truncated.truncate(tr)

    # Check the 3 untouched metadata entries match
    assert concatenate_multi_files_ts.meta.metadata[
        2:] == truncated.meta.metadata[1:]
    # Now check the truncated (but not truncated out) meta entry
    assert concatenate_multi_files_ts.meta.metadata[1][
        0].start != truncated.meta.metadata[0][0].start
    assert concatenate_multi_files_ts.meta.metadata[1][
        0].end == truncated.meta.metadata[0][0].end
    assert concatenate_multi_files_ts.meta.metadata[1][
        1] == truncated.meta.metadata[0][1]
    assert concatenate_multi_files_ts.meta.metadata[1][
        2] == truncated.meta.metadata[0][2]


def test_truncated_end_ts(concatenate_multi_files_ts):
    a = concatenate_multi_files_ts.meta.metadata[0][0].start - TimeDelta(1*u.day)
    b = concatenate_multi_files_ts.meta.metadata[-2][0].center
    tr = TimeRange(a, b)
    truncated = copy.deepcopy(concatenate_multi_files_ts)
    truncated = truncated.truncate(tr)

    # Check the 3 untouched metadata entries match
    assert concatenate_multi_files_ts.meta.metadata[:-2] == truncated.meta.metadata[:3]
    # Now check the truncated (but not truncated out) meta entry
    assert concatenate_multi_files_ts.meta.metadata[-2][
        0].start == truncated.meta.metadata[-1][0].start
    assert concatenate_multi_files_ts.meta.metadata[-2][
        0].end != truncated.meta.metadata[-1][0].end
    assert concatenate_multi_files_ts.meta.metadata[-2][
        1] == truncated.meta.metadata[-1][1]
    assert concatenate_multi_files_ts.meta.metadata[-2][
        2] == truncated.meta.metadata[-1][2]


def test_truncated_both_ts(concatenate_multi_files_ts):
    # This time range starts after and ends before the original, so expect truncation
    a = concatenate_multi_files_ts.meta.metadata[1][0].center
    b = concatenate_multi_files_ts.meta.metadata[-2][0].center
    tr = TimeRange(a, b)
    truncated = copy.deepcopy(concatenate_multi_files_ts)
    truncated = truncated.truncate(tr)

    # Check the 1 untouched metadata entry matches the original
    assert concatenate_multi_files_ts.meta.metadata[
        2:-2] == truncated.meta.metadata[1:-1]
    # Check the start truncated (but not truncated out) meta entry
    assert concatenate_multi_files_ts.meta.metadata[1][
        0].start != truncated.meta.metadata[0][0].start
    assert concatenate_multi_files_ts.meta.metadata[1][
        0].end == truncated.meta.metadata[0][0].end
    assert concatenate_multi_files_ts.meta.metadata[1][
        1] == truncated.meta.metadata[0][1]
    assert concatenate_multi_files_ts.meta.metadata[1][
        2] == truncated.meta.metadata[0][2]
    # Check the end truncated (but not truncated out) meta entry
    assert concatenate_multi_files_ts.meta.metadata[-2][
        0].start == truncated.meta.metadata[-1][0].start
    assert concatenate_multi_files_ts.meta.metadata[-2][
        0].end != truncated.meta.metadata[-1][0].end
    assert concatenate_multi_files_ts.meta.metadata[-2][
        1] == truncated.meta.metadata[-1][1]
    assert concatenate_multi_files_ts.meta.metadata[-2][
        2] == truncated.meta.metadata[-1][2]


@pytest.fixture
def truncated_new_tr_all_before_ts(concatenate_multi_files_ts):
    # Time range begins and ends before the data
    a = concatenate_multi_files_ts.meta.metadata[0][0].start - TimeDelta(2*u.day)
    b = concatenate_multi_files_ts.meta.metadata[0][0].start - TimeDelta(1*u.day)
    tr = TimeRange(a, b)
    truncated = copy.deepcopy(concatenate_multi_files_ts)
    truncated = truncated.truncate(tr)
    return truncated


@pytest.fixture
def truncated_new_tr_all_after_ts(concatenate_multi_files_ts):
    # Time range begins and ends after the data
    a = concatenate_multi_files_ts.meta.metadata[-1][0].end + TimeDelta(1*u.day)
    b = concatenate_multi_files_ts.meta.metadata[-1][0].end + TimeDelta(2*u.day)
    tr = TimeRange(a, b)
    truncated = copy.deepcopy(concatenate_multi_files_ts)
    truncated = truncated.truncate(tr)
    return truncated


def test_truncated_outside_tr_ts(truncated_new_tr_all_before_ts,
                                 truncated_new_tr_all_after_ts):
    assert (truncated_new_tr_all_before_ts.meta.metadata ==
            truncated_new_tr_all_after_ts.meta.metadata ==
            [])


def test_extraction(eve_test_ts):
    cmlon = eve_test_ts.extract('CMLon')
    # Test there's only one column in the data, metadata and units
    assert len(cmlon.to_dataframe().columns) == 1
    assert len(cmlon.meta.columns) == 1
    assert len(cmlon.units) == 1

    # Test this column name matches
    assert eve_test_ts.to_dataframe().columns[0] == eve_test_ts.to_dataframe().columns[0] == list(
        eve_test_ts.units.keys())[0]

    # Test the data matches
    extracted_df = DataFrame(eve_test_ts.to_dataframe()['CMLon']).dropna()
    extracted_df = extracted_df.sort_index()
    assert_frame_equal(cmlon.to_dataframe(), extracted_df)

# =============================================================================
# Test Concatenation Operations
# =============================================================================


@pytest.fixture
def concatenated_slices_test_ts(truncation_slice_test_ts_1, truncation_slice_test_ts_2,
                                truncation_slice_test_ts_3, truncation_slice_test_ts_4):
    # Concatenate the slices individually to make a TS similar to the original
    truncation_slice_test_ts_1 = truncation_slice_test_ts_1.concatenate(truncation_slice_test_ts_2)
    truncation_slice_test_ts_1 = truncation_slice_test_ts_1.concatenate(truncation_slice_test_ts_3)
    return truncation_slice_test_ts_1.concatenate(truncation_slice_test_ts_4)


@pytest.fixture
def concatenated_slices_test_list(truncation_slice_test_ts_1, truncation_slice_test_ts_2,
                                  truncation_slice_test_ts_3, truncation_slice_test_ts_4):
    # Concatenate the slices in a list to make a TS similar to the original
    return truncation_slice_test_ts_1.concatenate(
        [truncation_slice_test_ts_2, truncation_slice_test_ts_3, truncation_slice_test_ts_4]
    )


def test_concatenation_of_slices_ts(eve_test_ts, concatenated_slices_test_ts):
    # Test resulting DataFrame is similar to the original
    assert_frame_equal(concatenated_slices_test_ts.to_dataframe(), eve_test_ts.to_dataframe())
    # Otherwise: concatenated_ts.data.equals(eve_test_ts)
    # Compare timeranges from before and after match for both metadata and TS
    assert eve_test_ts.meta.time_range == concatenated_slices_test_ts.meta.time_range
    assert eve_test_ts.time_range == concatenated_slices_test_ts.time_range
    # Test metadata MetaDict matches
    eve_test_ts.meta.metadata[0][
        2] == concatenated_slices_test_ts.meta.metadata[0][
            2] == concatenated_slices_test_ts.meta.metadata[1][2]
    # ToDo: Will TSMD.concatenate() want to re-merge the metadata entries back into one?


def test_concatenation_of_slices_list(eve_test_ts, concatenated_slices_test_list):
    # Test resulting DataFrame is similar to the original
    assert_frame_equal(concatenated_slices_test_list.to_dataframe(), eve_test_ts.to_dataframe())
    # Otherwise: concatenated_ts.data.equals(eve_test_ts)
    # Compare timeranges from before and after match for both metadata and TS
    assert eve_test_ts.meta.time_range == concatenated_slices_test_list.meta.time_range
    assert eve_test_ts.time_range == concatenated_slices_test_list.time_range
    # Test metadata MetaDict matches
    eve_test_ts.meta.metadata[0][
        2] == concatenated_slices_test_list.meta.metadata[0][
            2] == concatenated_slices_test_list.meta.metadata[1][2]
    # ToDo: Will TSMD.concatenate() want to re-merge the metadata entries back into one?


@pytest.fixture
def different_data_concat(eve_test_ts, fermi_gbm_test_ts):
    # Take two different data sources and concatenate
    return eve_test_ts.concatenate(fermi_gbm_test_ts)


def test_concat_list(eve_test_ts, fermi_gbm_test_ts):
    assert (eve_test_ts.concatenate(fermi_gbm_test_ts) ==
            eve_test_ts.concatenate([fermi_gbm_test_ts]))


def test_concatenation_of_different_data_ts(eve_test_ts, fermi_gbm_test_ts,
                                            different_data_concat):
    # TODO: test the metadata is as expected using the below. (note ATM this fails if the order is changed)
    # assert different_data_concat.meta.metadata[0] == fermi_gbm_test_ts.meta.metadata[0]
    # assert different_data_concat.meta.metadata[1] == eve_test_ts.meta.metadata[0]
    value = True
    for key in list(different_data_concat.meta.metadata[0][2]
                    .keys()):
        if different_data_concat.meta.metadata[0][2][
                key] != fermi_gbm_test_ts.meta.metadata[0][2][key]:
            value = False
    for key in list(different_data_concat.meta.metadata[1][2]
                    .keys()):
        if different_data_concat.meta.metadata[1][2][
                key] != eve_test_ts.meta.metadata[0][2][key]:
            value = False
    assert value

    # Test units concatenation
    comined_units = copy.deepcopy(eve_test_ts.units)
    comined_units.update(fermi_gbm_test_ts.units)
    assert dict(different_data_concat.units) == dict(
        comined_units)

    # Test data is the concatenation
    comined_df = pd.concat([eve_test_ts.to_dataframe(), fermi_gbm_test_ts.to_dataframe()],
                           sort=False)
    comined_df = comined_df.sort_index()
    assert_frame_equal(different_data_concat.to_dataframe(), comined_df)


def test_concatenation_of_self(eve_test_ts):
    # Check that a self concatenation returns the original timeseries
    assert eve_test_ts.concatenate(eve_test_ts) == eve_test_ts
    assert eve_test_ts.concatenate([eve_test_ts]) == eve_test_ts


def test_concatenation_different_data_ts_error(eve_test_ts, fermi_gbm_test_ts):
    # Take two different data sources and concatenate but set with the same_source
    # kwarg as true. This should not concatenate.
    with pytest.raises(TypeError, match="TimeSeries classes must match if "
                                        "'same_source' is specified."):
        eve_test_ts.concatenate(fermi_gbm_test_ts, same_source=True)


def test_generic_construction_concatenation():
    nrows = 10
    # Generate the data and the corrisponding dates
    base = parse_time(datetime.datetime.today())
    times = base - TimeDelta(np.arange(nrows)*u.minute)
    intensity1 = np.sin(np.arange(0, 12 * np.pi, ((12 * np.pi) / (nrows))))
    intensity2 = np.sin(np.arange(0, 12 * np.pi, ((12 * np.pi) / (nrows))))

    # Create the data DataFrame, header MetaDict and units OrderedDict
    data = DataFrame(intensity1, index=times, columns=['intensity'])
    data2 = DataFrame(intensity2, index=times, columns=['intensity2'])
    units = OrderedDict([('intensity', u.W / u.m**2)])
    units2 = OrderedDict([('intensity2', u.W / u.m**2)])
    meta = MetaDict({'key': 'value'})
    meta2 = MetaDict({'key2': 'value2'})

    # Create TS individually
    ts_1 = sunpy.timeseries.TimeSeries(data, meta, units)
    ts_2 = sunpy.timeseries.TimeSeries(data2, meta2, units2)
    ts_concat = ts_1.concatenate(ts_2, axis=1)
    assert isinstance(ts_concat,
                      sunpy.timeseries.timeseriesbase.GenericTimeSeries)
    assert len(ts_concat.to_dataframe()) == len(times)
    assert ts_concat.columns == ['intensity', 'intensity2']
    assert len(ts_concat.meta.metadata) == 2

# =============================================================================
# Test Data Manipulation
# =============================================================================


@pytest.fixture
def column_quantity(eve_test_ts):
    # Get the astropy Quantity of values from a column
    return eve_test_ts.quantity('CMLon')


def test_column_quantity(eve_test_ts, column_quantity):
    # Test the units and values match
    assert eve_test_ts.units['CMLon'] == column_quantity.unit
    assert ((eve_test_ts.to_dataframe()['CMLon'].values == column_quantity.value) |
            (np.isnan(eve_test_ts.to_dataframe()['CMLon'].values) &
             np.isnan(column_quantity.value))).all()


def test_add_column_from_quantity(eve_test_ts, column_quantity):
    new_ts = eve_test_ts.add_column('quantity_added', column_quantity)
    # Test the column similar to the original quantity?
    assert_quantity_allclose(new_ts.quantity('quantity_added'), column_quantity)
    # Test the full list of columns are pressent
    assert set(new_ts.to_dataframe().columns) == set(
        eve_test_ts.to_dataframe().columns) | {'quantity_added'}


def test_remove_column(eve_test_ts):
    removed = eve_test_ts.remove_column('XRS-B proxy')
    # Check that column remains in eve_test_ts but isn't in removed
    assert 'XRS-B proxy' not in removed.columns
    assert 'XRS-B proxy' in eve_test_ts.columns
    assert len(removed.columns) == len(eve_test_ts.columns) - 1 == 18

    # Check units updated correctly
    assert len(removed.columns) == len(removed.units)

    # Check data updated correctly
    assert len(removed.columns) == removed.to_dataframe().shape[1]

    # Check that removing a non-existant column errors
    with pytest.raises(ValueError):
        eve_test_ts.remove_column('random column name')


def test_add_column_from_array(eve_test_ts, column_quantity):
    # Test the column similar to the original quantity?
    new_ts = eve_test_ts.add_column(
        'array_added', column_quantity.value, unit=column_quantity.unit)
    assert_quantity_allclose(
        new_ts.quantity('array_added'), column_quantity)

    # Test the full list of columns are pressent
    assert set(new_ts.to_dataframe().columns) == set(
        eve_test_ts.to_dataframe().columns) | {'array_added'}


def test_add_column_from_array_no_units(eve_test_ts, column_quantity):
    ts = eve_test_ts.add_column('array_added', column_quantity.value)
    assert (ts.quantity('array_added') == column_quantity.value).all()

# =============================================================================
# Test Exporting to different formats
# =============================================================================


def test_ts_to_table(generic_ts):
    tbl = generic_ts.to_table()
    assert isinstance(tbl, Table)
    assert tbl.keys() == ['date', *generic_ts.columns]
    assert len(tbl) == len(generic_ts.to_dataframe())
    assert (tbl[generic_ts.columns[0]].quantity ==
            generic_ts.quantity(generic_ts.columns[0])).all()


def test_ts_to_dataframe(generic_ts):
    df = generic_ts.to_dataframe()
    assert isinstance(df, DataFrame)
    assert_frame_equal(df, generic_ts.to_dataframe())


def test_ts_to_array(generic_ts):
    arr = generic_ts.to_array()
    assert isinstance(arr, np.ndarray)
    assert len(arr) == len(generic_ts.to_dataframe())

# =============================================================================
# Test Basic Working Peek
# =============================================================================


@figure_test
def test_eve_peek(eve_test_ts):
    eve_test_ts.peek()


@figure_test
def test_esp_peek(esp_test_ts):
    esp_test_ts.peek()


# This warning is fixed in matplotlib, and the filter can be removed once
# matplotlib 3.3.1 is released (https://github.com/matplotlib/matplotlib/pull/18101)
@pytest.mark.filterwarnings('ignore:Support for multi-dimensional indexing.*is deprecated')
@figure_test
def test_fermi_gbm_peek(fermi_gbm_test_ts):
    fermi_gbm_test_ts.peek()


@figure_test
def test_norh_peek(norh_test_ts):
    norh_test_ts.peek()


@figure_test
def test_goes_peek(goes_test_ts):
    goes_test_ts.peek()


@figure_test
def test_lyra_peek(lyra_test_ts):
    lyra_test_ts.peek()


@figure_test
def test_rhessi_peek(rhessi_test_ts):
    rhessi_test_ts.peek()


@figure_test
def test_noaa_json_ind_peek(noaa_ind_json_test_ts):
    noaa_ind_json_test_ts.peek()


@figure_test
def test_noaa_json_pre_peek(noaa_pre_json_test_ts):
    noaa_pre_json_test_ts.peek()


@figure_test
def test_generic_ts_peek(generic_ts):
    generic_ts.peek()


@figure_test
def test_column_subset_peek(generic_ts):
    generic_ts.peek(columns=['intensity2'])


def test_empty_ts_invalid_peek(generic_ts):
    # Truncate a timeseries so it's empty
    a = generic_ts.time_range.start - TimeDelta(2*u.day)
    b = generic_ts.time_range.start - TimeDelta(1*u.day)
    empty_ts = generic_ts.truncate(TimeRange(a, b))
    with pytest.raises(ValueError):
        empty_ts.peek()

# =============================================================================
# Test Other Functions
# =============================================================================


def test_equality(generic_ts, table_ts):
    generic_copy_ts = copy.deepcopy(generic_ts)
    assert generic_ts == generic_copy_ts
    generic_copy_ts.meta.metadata[0][2]['key'] = 1
    assert generic_ts != generic_copy_ts
    assert generic_ts != table_ts


def test_equality_different_ts_types(generic_ts):
    # this should fail as they're not the smae type and can't match
    assert not (generic_ts == eve_test_ts)


def test_ts_index(generic_ts):
    assert (generic_ts.index == generic_ts.to_dataframe().index).all()


def test_ts_shape(generic_ts):
    assert generic_ts.shape == generic_ts.to_dataframe().shape


def test_ts_sort_index(generic_ts):
    assert generic_ts.sort_index().to_dataframe().equals(generic_ts.to_dataframe().sort_index())


def test_timeseries_array():
    times = parse_time("now") - TimeDelta(np.arange(5) * u.minute)
    intensity = np.sin(np.arange(0, 12 * np.pi, step=(12 * np.pi) / (5)))
    data = np.array([times, intensity]).T
    with pytest.warns(SunpyUserWarning, match='Unknown units'):
        ts = sunpy.timeseries.TimeSeries(data, {})
        assert isinstance(ts, sunpy.timeseries.GenericTimeSeries)


# TODO:
# _validate_units
# _validate_meta
# Extracting column as quantity or array#ts_eve = ts_eve.add_column(colname, qua_new, overwrite=True)
# Updating a column using quantity or array#ts_eve = ts_eve.add_column(colname, qua_new, overwrite=True)
# Updating the units# ts_eve = ts_eve.add_column(colname, qua_new, unit=unit, overwrite=True)
