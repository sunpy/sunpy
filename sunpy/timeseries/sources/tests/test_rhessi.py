import sunpy.timeseries
from sunpy.data.test import get_test_filepath

rhessi_filepath = get_test_filepath('hsi_obssumm_20120601_018_truncated.fits.gz')


def test_implicit_rhessi():
    # Test a RHESSI TimeSeries
    ts_rhessi = sunpy.timeseries.TimeSeries(rhessi_filepath)
    assert isinstance(ts_rhessi, sunpy.timeseries.sources.rhessi.RHESSISummaryTimeSeries)


def test_rhessi():
    # Test a RHESSI TimeSeries
    ts_rhessi = sunpy.timeseries.TimeSeries(rhessi_filepath, source='RHESSI')
    assert isinstance(ts_rhessi, sunpy.timeseries.sources.rhessi.RHESSISummaryTimeSeries)
