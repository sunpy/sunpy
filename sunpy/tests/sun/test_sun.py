from __future__ import absolute_import

from datetime import datetime

import sunpy.sun as sun
from numpy.testing import assert_array_almost_equal

def test_sunearth_distance():
    assert_array_almost_equal(sun.sunearth_distance("2010/02/04"), 0.9858, decimal=4)
    assert_array_almost_equal(sun.sunearth_distance("2009/04/13"), 1.003, decimal=4)
    assert_array_almost_equal(sun.sunearth_distance("2008/06/20"), 1.016, decimal=4)
    assert_array_almost_equal(sun.sunearth_distance("2007/08/15"), 1.013, decimal=4)
    assert_array_almost_equal(sun.sunearth_distance("2007/10/02"), 1.001, decimal=4)
    assert_array_almost_equal(sun.sunearth_distance("2006/12/27"), 0.9834, decimal=4)


