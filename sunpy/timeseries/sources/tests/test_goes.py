import pytest

import sunpy.timeseries
from sunpy.data.test import get_test_filepath
from sunpy.tests.helpers import figure_test
from sunpy.util.exceptions import SunpyUserWarning

goes_filepath_com = get_test_filepath('go1520120601.fits.gz')
new_goes15_filepath = get_test_filepath('goes_truncated_test_goes15.nc')
new_goes17_filepath = get_test_filepath('goes_truncated_test_goes17.nc')
goes13_leap_second_filepath = get_test_filepath('goes_13_leap_second.nc')


def test_implicit_goes(goes_test_ts):
    assert isinstance(goes_test_ts, sunpy.timeseries.sources.goes.XRSTimeSeries)


def test_implicit_goes_com(goes_test_ts):
    assert isinstance(goes_test_ts, sunpy.timeseries.sources.goes.XRSTimeSeries)


def test_implicit_new_goes15():
    # Test a GOES TimeSeries
    ts_goes = sunpy.timeseries.TimeSeries(new_goes15_filepath)
    assert isinstance(ts_goes, sunpy.timeseries.sources.goes.XRSTimeSeries)


def test_implicit_new_goes17():
    # Test a GOES TimeSeries
    ts_goes = sunpy.timeseries.TimeSeries(new_goes17_filepath)
    assert isinstance(ts_goes, sunpy.timeseries.sources.goes.XRSTimeSeries)


def test_implicit_goes_satno(goes_test_ts):
    assert goes_test_ts.observatory == 'GOES-15'


def test_implicit_new_goes15_satno():
    # Test a GOES TimeSeries for satellite number
    ts_goes = sunpy.timeseries.TimeSeries(new_goes15_filepath)
    assert ts_goes.observatory == 'GOES-15'


def test_implicit_new_goes17_satno():
    # Test a GOES TimeSeries for satellite number
    ts_goes = sunpy.timeseries.TimeSeries(new_goes17_filepath)
    assert ts_goes.observatory == 'GOES-17'


def test_implicit_goes_satno_missing():
    # Test a GOES TimeSeries for a missing satellite number
    ts_goes = sunpy.timeseries.TimeSeries(new_goes17_filepath)
    del ts_goes.meta.metas[0]['id']
    assert ts_goes.observatory is None


def test_goes_com():
    # Test a GOES TimeSeries
    ts_goes = sunpy.timeseries.TimeSeries(goes_filepath_com, source='XRS')
    assert isinstance(ts_goes, sunpy.timeseries.sources.goes.XRSTimeSeries)


def test_new_goes15():
    # Test a GOES TimeSeries
    ts_goes = sunpy.timeseries.TimeSeries(new_goes15_filepath, source='XRS')
    assert isinstance(ts_goes, sunpy.timeseries.sources.goes.XRSTimeSeries)


def test_new_goes16():
    # Test a GOES TimeSeries
    ts_goes = sunpy.timeseries.TimeSeries(new_goes17_filepath, source='XRS')
    assert isinstance(ts_goes, sunpy.timeseries.sources.goes.XRSTimeSeries)


def test_goes_netcdf_time_parsing15():
    # testing to make sure the time is correctly parsed
    ts_goes = sunpy.timeseries.TimeSeries(new_goes15_filepath, source="XRS")
    assert ts_goes.time[0].isot == '2013-10-28T00:00:01.385'


def test_goes_netcdf_time_parsing17():
    # testing to make sure the time is correctly parsed
    ts_goes = sunpy.timeseries.TimeSeries(new_goes17_filepath, source="XRS")
    assert ts_goes.time[0].isot == '2020-10-16T00:00:00.477'


@pytest.mark.remote_data
def test_goes_remote():
    # Older format file
    goes = sunpy.timeseries.TimeSeries(
        'https://umbra.nascom.nasa.gov/goes/fits/1986/go06860129.fits')
    assert isinstance(goes, sunpy.timeseries.sources.goes.XRSTimeSeries)
    # Newer format
    goes = sunpy.timeseries.TimeSeries(
        'https://umbra.nascom.nasa.gov/goes/fits/2018/go1520180626.fits')
    assert isinstance(goes, sunpy.timeseries.sources.goes.XRSTimeSeries)


def test_goes_leap_seconds():
    with pytest.warns(SunpyUserWarning, match="There is one leap second timestamp present in: goes_13_leap_second"):
        ts = sunpy.timeseries.TimeSeries(goes13_leap_second_filepath)
    assert ts.time[-1].isot == '2015-06-30T23:59:59.999'


def test_goes_plot_column(goes_test_ts):
    ax = goes_test_ts.plot(columns=['xrsa'])
    assert len(ax.lines) == 1
    assert '0.5--4.0' == ax.lines[0].get_label().split()[0]


@figure_test
def test_goes_peek(goes_test_ts):
    goes_test_ts.peek()
