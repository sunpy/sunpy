import sunpy.timeseries
from sunpy.data.test import get_test_filepath

norh_filepath = get_test_filepath('tca110810_truncated')


def test_implicit_norh():
    # Test a NoRH TimeSeries
    ts_norh = sunpy.timeseries.TimeSeries(norh_filepath)
    assert isinstance(ts_norh, sunpy.timeseries.sources.norh.NoRHTimeSeries)


def test_norh():
    # Test a NoRH TimeSeries
    ts_norh = sunpy.timeseries.TimeSeries(norh_filepath, source='NoRH')
    assert isinstance(ts_norh, sunpy.timeseries.sources.norh.NoRHTimeSeries)
