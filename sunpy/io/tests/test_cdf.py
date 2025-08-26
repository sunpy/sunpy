
import numpy as np
import pytest

import astropy.units as u

from sunpy.data.test import get_test_filepath
from sunpy.io._cdf import read_cdf
from sunpy.net import Fido
from sunpy.net import attrs as a
from sunpy.timeseries import GenericTimeSeries, TimeSeries
from sunpy.util.exceptions import SunpyUserWarning

filepath = get_test_filepath('solo_L2_epd-ept-north-hcad_20200713_V02.cdf')


def test_read_cdf():
    all_ts = read_cdf(filepath)
    assert isinstance(all_ts, list)
    assert len(all_ts) == 3

    ts = all_ts[0]
    assert isinstance(ts, GenericTimeSeries)

    col = ts.quantity('Electron_Flux_0')
    assert col.unit == u.Unit("1 / (cm2 MeV s sr)")
    # Check that fillvals are replaced by NaN
    assert np.sum(np.isnan(col)) == 189


@pytest.mark.remote_data
def test_read_psp_data():
    # This was a failing example provided by
    # https://github.com/sunpy/sunpy/issues/7565
    dataset = 'PSP_SWP_SPI_SF00_L3_MOM'
    trange = a.Time('2023-03-14', '2023-03-15')
    result = Fido.search(trange, a.cdaweb.Dataset(dataset))
    downloaded_files = Fido.fetch(result)

    with pytest.warns(SunpyUserWarning, match="astropy did not recognize units of"):
        ts = TimeSeries(downloaded_files, concatenate=True)

    assert isinstance(ts, GenericTimeSeries)
    col = ts.quantity('EFLUX_VS_ENERGY_0')
    assert np.sum(np.isnan(col)) >= 0
