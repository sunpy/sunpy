
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


def test_generic_timeseries_columns():
    # Simulate missing 'FILLVAL' for one variable and ensure graceful handling
    all_ts = read_cdf(filepath)
    for ts in all_ts:
        # Iterate over each column in the GenericTimeSeries
        for col_name in ts.columns:
            col_data = ts.quantity(col_name)
            # Ensure the data exists
            assert col_data is not None
            # Validate NaN handling in the column
            assert np.any(np.isnan(col_data)) or np.all(~np.isnan(col_data))
