
import numpy as np

import astropy.units as u

from sunpy.data.test import get_test_filepath
from sunpy.io._cdf import read_cdf
from sunpy.timeseries import GenericTimeSeries

filepath = get_test_filepath('solo_L2_epd-ept-north-hcad_20200713_V02.cdf')


def test_read_cdf():
    all_ts = read_cdf(filepath)
    assert isinstance(all_ts, list)
    assert len(all_ts) == 3

    ts = all_ts[0]
    print(ts.columns)
    assert isinstance(ts, GenericTimeSeries)

    col = ts.quantity('Electron_Flux_0')
    assert col.unit == u.Unit("1 / (cm2 MeV s sr)")
    # Check that fillvals are replaced by NaN
    assert np.sum(np.isnan(col)) == 189

def test_check_nan_values():
    all_ts = read_cdf(filepath)
    assert isinstance(all_ts, list)

    assert len(all_ts) > 0

    for i, ts in enumerate(all_ts):
        print(f"Checking Timeseries {i + 1} for NaN values...")
        has_nan = False

        data = ts.to_dataframe()
        for column in data.columns:
            if data[column].isna().any():
                has_nan = True
                print(f"  Column '{column}' contains NaN values.")

        if not has_nan:
            print(f"  No NaN values found in Timeseries {i + 1}.")
        print()
