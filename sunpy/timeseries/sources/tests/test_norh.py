import sunpy.timeseries
from sunpy.data.test import get_test_filepath
from sunpy.tests.helpers import figure_test

norh_filepath = get_test_filepath('tca110810_truncated')


def test_implicit_norh():
    # Test a NoRH TimeSeries
    ts_norh = sunpy.timeseries.TimeSeries(norh_filepath)
    assert isinstance(ts_norh, sunpy.timeseries.sources.norh.NoRHTimeSeries)


def test_norh():
    # Test a NoRH TimeSeries
    ts_norh = sunpy.timeseries.TimeSeries(norh_filepath, source='NoRH')
    assert isinstance(ts_norh, sunpy.timeseries.sources.norh.NoRHTimeSeries)


def test_norh_plot_column(norh_test_ts):
    ax = norh_test_ts.plot(columns=['Correlation Coefficient'])
    assert len(ax.lines) == 1
    assert '17GHZ' == ax.lines[0].get_label()


@figure_test
def test_norh_peek(norh_test_ts):
    norh_test_ts.peek()
