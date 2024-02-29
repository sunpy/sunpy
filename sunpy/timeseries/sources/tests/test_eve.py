
import pytest

import sunpy.timeseries
from sunpy.data.test import get_test_filepath
from sunpy.tests.helpers import figure_test
from sunpy.util.exceptions import SunpyUserWarning

esp_filepath = get_test_filepath('eve_l1_esp_2011046_00_truncated.fits')
eve_filepath = get_test_filepath('EVE_L0CS_DIODES_1m_truncated.txt')


@pytest.mark.skipif(pytest.__version__ < "8.0.0", reason="pytest >= 8.0.0 raises an extra warning for this test")
def test_eve():
    with pytest.warns(SunpyUserWarning, match='Unknown units for oldXRSB proxy'):
        with pytest.warns(SunpyUserWarning, match='Unknown units for x_cool'):
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


def test_esp_plot_column(esp_test_ts):
    axes = esp_test_ts.plot(columns=['QD', 'CH_18', 'CH_36'])
    assert len(axes) == 3
    assert '0.1-7nm' in axes[0].get_ylabel()
    assert '18nm' in axes[1].get_ylabel()
    assert '36nm' in axes[2].get_ylabel()


@figure_test
def test_esp_peek(esp_test_ts):
    esp_test_ts.peek()


@figure_test
def test_eve_peek(eve_test_ts):
    eve_test_ts.peek()
