from pathlib import Path

import numpy as np

from sunpy.data.test import get_dummy_map_from_header, get_test_filepath

TEST_IRIS_SJI_HEADER = get_test_filepath('iris_l2_20130801_074720_4040000014_SJI_1400_t000.header')

def test_get_dummy_map_from_header_b_scale_zero_applied():
    # We use the IRIS one since it contains a BZERO and BSCALE in the test header
    iris_map = get_dummy_map_from_header(Path(TEST_IRIS_SJI_HEADER))
    assert iris_map.meta["bscale"] == 0.25
    assert iris_map.meta["bzero"] == 7992
    # Data from the test header is scaled by BSCALE and BZERO
    # but the max value it will return is 100 and converted to int
    # So its 24 instead of 25 due to the rounding
    assert np.max(iris_map.data) == 8016
