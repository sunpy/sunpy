import pytest

import sunpy.timeseries
from sunpy.data.test import get_test_filepath

goes_filepath = get_test_filepath('go1520110607.fits')
goes_filepath_com = get_test_filepath('go1520120601.fits.gz')
new_goes15_filepath = get_test_filepath('goes_truncated_test_goes15.nc')
new_goes17_filepath = get_test_filepath('goes_truncated_test_goes17.nc')


def test_implicit_goes():
    # Test a GOES TimeSeries
    ts_goes = sunpy.timeseries.TimeSeries(goes_filepath)
    assert isinstance(ts_goes, sunpy.timeseries.sources.goes.XRSTimeSeries)


def test_implicit_goes_com():
    # Test a GOES TimeSeries
    ts_goes = sunpy.timeseries.TimeSeries(goes_filepath_com)
    assert isinstance(ts_goes, sunpy.timeseries.sources.goes.XRSTimeSeries)


def test_implicit_new_goes15():
    # Test a GOES TimeSeries
    ts_goes = sunpy.timeseries.TimeSeries(new_goes15_filepath)
    assert isinstance(ts_goes, sunpy.timeseries.sources.goes.XRSTimeSeries)


def test_implicit_new_goes17():
    # Test a GOES TimeSeries
    ts_goes = sunpy.timeseries.TimeSeries(new_goes17_filepath)
    assert isinstance(ts_goes, sunpy.timeseries.sources.goes.XRSTimeSeries)


def test_implicit_goes_satno():
    # Test a GOES TimeSeries for satellite number
    ts_goes = sunpy.timeseries.TimeSeries(goes_filepath)
    assert ts_goes.observatory == 'GOES-15'


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


def test_goes():
    # Test a GOES TimeSeries
    ts_goes = sunpy.timeseries.TimeSeries(goes_filepath, source='XRS')
    assert isinstance(ts_goes, sunpy.timeseries.sources.goes.XRSTimeSeries)


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
    # testing to make sure the time is correctly parsed (to ignore leap seconds)
    ts_goes = sunpy.timeseries.TimeSeries(new_goes15_filepath, source="XRS")
    assert ts_goes.index[0].strftime("%Y-%m-%d %H:%M:%S.%f") == '2013-10-28 00:00:01.385000'


def test_goes_netcdf_time_parsing17():
    # testing to make sure the time is correctly parsed (to ignore leap seconds)
    ts_goes = sunpy.timeseries.TimeSeries(new_goes17_filepath, source="XRS")
    assert ts_goes.index[0].strftime("%Y-%m-%d %H:%M:%S.%f") == '2020-10-16 00:00:00.476771'


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
