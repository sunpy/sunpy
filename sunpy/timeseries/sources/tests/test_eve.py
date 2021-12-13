import pytest

import sunpy.timeseries
from sunpy.data.test import get_test_filepath
from sunpy.util.exceptions import SunpyUserWarning

esp_filepath = get_test_filepath('eve_l1_esp_2011046_00_truncated.fits')
eve_filepath = get_test_filepath('EVE_L0CS_DIODES_1m_truncated.txt')


def test_eve():
    # Test an EVE TimeSeries
    with pytest.warns(SunpyUserWarning, match='Unknown units for x_cool proxy'):
        ts_eve = sunpy.timeseries.TimeSeries(eve_filepath, source='EVE')
    assert isinstance(ts_eve, sunpy.timeseries.sources.eve.EVESpWxTimeSeries)


def test_implicit_esp():
    # Test an ESP TimeSeries
    ts_esp = sunpy.timeseries.TimeSeries(esp_filepath)
    assert isinstance(ts_esp, sunpy.timeseries.sources.eve.ESPTimeSeries)


def test_esp():
    # Test an ESP TimeSeries
    ts_esp = sunpy.timeseries.TimeSeries(esp_filepath, source='ESP')
    assert isinstance(ts_esp, sunpy.timeseries.sources.eve.ESPTimeSeries)
