import sunpy.timeseries
from sunpy.data.test import get_test_filepath
from sunpy.tests.helpers import figure_test

rhessi_filepath = get_test_filepath('hsi_obssumm_20120601_018_truncated.fits.gz')


def test_implicit_rhessi():
    # Test a RHESSI TimeSeries
    ts_rhessi = sunpy.timeseries.TimeSeries(rhessi_filepath)
    assert isinstance(ts_rhessi, sunpy.timeseries.sources.rhessi.RHESSISummaryTimeSeries)


def test_rhessi():
    # Test a RHESSI TimeSeries
    ts_rhessi = sunpy.timeseries.TimeSeries(rhessi_filepath, source='RHESSI')
    assert isinstance(ts_rhessi, sunpy.timeseries.sources.rhessi.RHESSISummaryTimeSeries)


def test_rhessi_plot_column(rhessi_test_ts):
    ax = rhessi_test_ts.plot()
    assert len(ax.lines) == 9
    assert '3 - 6 keV' == ax.lines[0].get_label()
    assert '6 - 12 keV' == ax.lines[1].get_label()
    assert '12 - 25 keV' == ax.lines[2].get_label()


@figure_test
def test_rhessi_peek(rhessi_test_ts):
    rhessi_test_ts.peek()
