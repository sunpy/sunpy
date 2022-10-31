
import sunpy.timeseries
from sunpy.data.test import get_test_filepath
from sunpy.tests.helpers import figure_test

fermi_gbm_filepath = get_test_filepath('gbm.fits')


def test_implicit_fermi_gbm():
    # Test a GBMSummary TimeSeries
    ts_gbm = sunpy.timeseries.TimeSeries(fermi_gbm_filepath)
    assert isinstance(ts_gbm, sunpy.timeseries.sources.fermi_gbm.GBMSummaryTimeSeries)


def test_fermi_gbm():
    # Test a GBMSummary TimeSeries
    ts_gbm = sunpy.timeseries.TimeSeries(fermi_gbm_filepath, source='GBMSummary')
    assert isinstance(ts_gbm, sunpy.timeseries.sources.fermi_gbm.GBMSummaryTimeSeries)


def test_fermi_gbm_plot_column(fermi_gbm_test_ts):
    ax = fermi_gbm_test_ts.plot(columns=['4-15 keV', '100-300 keV'])
    assert len(ax.lines) == 2
    assert '4-15 keV' == ax.lines[0].get_label()
    assert '100-300 keV' == ax.lines[1].get_label()


@figure_test
def test_fermi_gbm_peek(fermi_gbm_test_ts):
    fermi_gbm_test_ts.peek()
