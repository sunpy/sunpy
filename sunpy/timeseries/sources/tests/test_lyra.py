import sunpy.timeseries
from sunpy.data.test import get_test_filepath

lyra_filepath = get_test_filepath('lyra_20150101-000000_lev3_std_truncated.fits.gz')


def test_implicit_lyra():
    # Test a LYRA TimeSeries
    ts_lyra = sunpy.timeseries.TimeSeries(lyra_filepath)
    assert isinstance(ts_lyra, sunpy.timeseries.sources.lyra.LYRATimeSeries)


def test_lyra():
    # Test a LYRA TimeSeries
    ts_lyra = sunpy.timeseries.TimeSeries(lyra_filepath, source='LYRA')
    assert isinstance(ts_lyra, sunpy.timeseries.sources.lyra.LYRATimeSeries)
