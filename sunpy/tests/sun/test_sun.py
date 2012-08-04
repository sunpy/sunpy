from __future__ import absolute_import

from sunpy.sun import sun
from numpy.testing import assert_array_almost_equal


def test_sunearth_distance():
    # Source for these values 
    # wolframalpha.com
    # http://www.wolframalpha.com/input/?i=earth-sun+distance+on+2010%2F02%2F04
    assert_array_almost_equal(sun.sunearth_distance("2010/02/04"), 0.9858, decimal=3)
    assert_array_almost_equal(sun.sunearth_distance("2009/04/13"), 1.003, decimal=3)
    assert_array_almost_equal(sun.sunearth_distance("2008/06/20"), 1.016, decimal=3)
    assert_array_almost_equal(sun.sunearth_distance("2007/08/15"), 1.013, decimal=3)
    assert_array_almost_equal(sun.sunearth_distance("2007/10/02"), 1.001, decimal=3)
    assert_array_almost_equal(sun.sunearth_distance("2006/12/27"), 0.9834, decimal=3)


