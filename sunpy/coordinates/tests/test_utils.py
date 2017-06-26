from __future__ import absolute_import, division, print_function

import pytest

import numpy as np
from astropy.coordinates import SkyCoord
import astropy.units as u
import sunpy.coordinates
import sunpy.map
import sunpy.data.test
from sunpy.coordinates.great_arc import great_arc, _calculate_great_arc


# Test the great arc code
def test_great_arc():
    # Number of points in the return
    num = 3

    m = sunpy.map.Map(sunpy.map.Map(sunpy.data.test.get_test_filepath('aia_171_level1.fits')))
    coordinate_frame = m.coordinate_frame
    a = SkyCoord(600*u.arcsec, -600*u.arcsec, frame=coordinate_frame)
    b = SkyCoord(-100*u.arcsec, 800*u.arcsec, frame=coordinate_frame)
    v = great_arc(a, b, number_points=num)

    # Return a SkyCoord
    assert isinstance(v, SkyCoord)

    # Correct number of points along the great arc
    assert len(v) == 3

    # Start and end point values are as expected
    np.testing.assert_almost_equal(v[0].Tx.value, a.Tx.value)
    np.testing.assert_almost_equal(v[0].Ty.value, a.Ty.value)
    np.testing.assert_almost_equal(v[-1].Tx.value, b.Tx.value)
    np.testing.assert_almost_equal(v[-1].Ty.value, b.Ty.value)

    # Make sure the output observer is correct
    for i in range(0, len(v)):
        assert v[i].observer == a.observer


# Test the calculation of the great arc.  Different centers - zero and non-zero.
# Tests that the great arc calculation correctly accounts for the location of
# the center.
@pytest.mark.parametrize("center", [np.asarray([0, 0, 0]), np.asarray([1, 2, 3])])
def test_calculate_great_arc(center):
    # Testing accuracy
    decimal = 6

    # Number of points in the return
    num = 3

    # Make sure everything works when z is zero.
    a = np.asarray([1, 0, 0])
    b = np.asarray([0, 1, 0])
    test_a = a + center
    test_b = b + center
    v_xyz = _calculate_great_arc(test_a, test_b, center, num)
    assert v_xyz.shape == (3, 3)
    np.testing.assert_almost_equal(v_xyz[0, :], test_a, decimal=decimal)
    np.testing.assert_almost_equal(v_xyz[1, :], np.asarray([7.07106781e-01, 7.07106781e-01, 0.0]) + center, decimal=decimal)
    np.testing.assert_almost_equal(v_xyz[2, :], test_b, decimal=decimal)

    # Make sure everything works when z is non-zero.
    c = np.asarray([1, 0, 1])
    d = np.asarray([0, 1, 1])
    test_c = c + center
    test_d = d + center
    v_xyz = _calculate_great_arc(test_c, test_d, center, num)
    np.testing.assert_almost_equal(v_xyz[0, :], test_c, decimal=decimal)
    np.testing.assert_almost_equal(v_xyz[1, :], np.asarray([5.77350269e-01, 5.77350269e-01, 1.15470054e+00]) + center, decimal=decimal)
    np.testing.assert_almost_equal(v_xyz[2, :], test_d, decimal=decimal)
