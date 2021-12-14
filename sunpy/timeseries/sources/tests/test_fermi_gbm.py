import sunpy.timeseries
from sunpy.data.test import get_test_filepath

fermi_gbm_filepath = get_test_filepath('gbm.fits')


def test_implicit_fermi_gbm():
    # Test a GBMSummary TimeSeries
    ts_gbm = sunpy.timeseries.TimeSeries(fermi_gbm_filepath)
    assert isinstance(ts_gbm, sunpy.timeseries.sources.fermi_gbm.GBMSummaryTimeSeries)


def test_fermi_gbm():
    # Test a GBMSummary TimeSeries
    ts_gbm = sunpy.timeseries.TimeSeries(fermi_gbm_filepath, source='GBMSummary')
    assert isinstance(ts_gbm, sunpy.timeseries.sources.fermi_gbm.GBMSummaryTimeSeries)
