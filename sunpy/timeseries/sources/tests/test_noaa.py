import pytest

import sunpy.timeseries
from sunpy.data.test import get_test_filepath
from sunpy.tests.helpers import figure_test

noaa_ind_json_filepath = get_test_filepath('observed-solar-cycle-indices-truncated.json')
noaa_pre_json_filepath = get_test_filepath('predicted-solar-cycle-truncated.json')


def test_noaa_ind_json():
    # Test a NOAAPredictIndices TimeSeries json
    ts_noaa_ind = sunpy.timeseries.TimeSeries(noaa_ind_json_filepath, source='NOAAIndices')
    assert isinstance(ts_noaa_ind, sunpy.timeseries.sources.noaa.NOAAIndicesTimeSeries)


# The pre- data involves dates long in the future, so ignore an ERFA warning
# when parsing these dates.
@pytest.mark.filterwarnings('ignore:ERFA function.*dubious year')
def test_noaa_pre_json():
    # Test a NOAAIndices TimeSeries json
    ts_noaa_pre = sunpy.timeseries.TimeSeries(
        noaa_pre_json_filepath, source='NOAAPredictIndices')
    assert isinstance(ts_noaa_pre, sunpy.timeseries.sources.noaa.NOAAPredictIndicesTimeSeries)


def test_noaa_json_pre_plot_column(noaa_pre_json_test_ts):
    ax = noaa_pre_json_test_ts.plot(columns=['sunspot', 'sunspot high', 'sunspot low'])
    assert len(ax.lines) == 3
    assert 'sunspot' == ax.lines[0].get_label()
    assert 'sunspot high' == ax.lines[1].get_label()
    assert 'sunspot low' == ax.lines[2].get_label()


def test_noaa_json_ind_plot_column(noaa_ind_json_test_ts):
    ax = noaa_ind_json_test_ts.plot(columns=['sunspot SWO', 'sunspot SWO smooth'])
    assert len(ax.lines) == 2
    assert 'sunspot SWO' == ax.lines[0].get_label()
    assert 'sunspot SWO smooth' == ax.lines[1].get_label()


@figure_test
def test_noaa_json_pre_peek(noaa_pre_json_test_ts):
    noaa_pre_json_test_ts.peek()


@figure_test
def test_noaa_json_ind_peek(noaa_ind_json_test_ts):
    noaa_ind_json_test_ts.peek()
