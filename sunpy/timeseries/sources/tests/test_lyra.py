import sunpy.timeseries
from sunpy.data.test import get_test_filepath
from sunpy.tests.helpers import figure_test

lyra_filepath = get_test_filepath('lyra_20150101-000000_lev3_std_truncated.fits.gz')


def test_implicit_lyra():
    # Test a LYRA TimeSeries
    ts_lyra = sunpy.timeseries.TimeSeries(lyra_filepath)
    assert isinstance(ts_lyra, sunpy.timeseries.sources.lyra.LYRATimeSeries)


def test_lyra():
    # Test a LYRA TimeSeries
    ts_lyra = sunpy.timeseries.TimeSeries(lyra_filepath, source='LYRA')
    assert isinstance(ts_lyra, sunpy.timeseries.sources.lyra.LYRATimeSeries)


def test_lyra_plot_column(lyra_test_ts):
    axes = lyra_test_ts.plot(columns=['CHANNEL1', 'CHANNEL3'])
    assert len(axes) == 2
    assert 'Lyman alpha' in axes[0].get_ylabel()
    assert 'Al filter' in axes[1].get_ylabel()


@figure_test
def test_lyra_peek(lyra_test_ts):
    lyra_test_ts.peek()
