from __future__ import absolute_import, division, print_function

import pytest

import numpy as np
from astropy.coordinates import SkyCoord
import astropy.units as u
from sunpy.coordinates import frames
import sunpy.map
import sunpy.data.test
from sunpy.coordinates.utils import GreatArc


# Test the great arc code
# Test the calculation of coordinates using varying numbers of points on
# initialization of the GreatArc object
@pytest.mark.parametrize("points_requested, points_expected",
                         [(None, 100),  # Test default
                          (3, 3),  # Test int as an option
                          (np.linspace(0, 1, 43), 43),  # Test equally spaced monotonically increasing numpy array
                          (np.asarray([0.1, 0.2, 0.6, 0.67, 0.99]), 5),  # Test unequally spaced monotonically increasing numpy array
                          (np.asarray([0.93, 0.78, 0.3, 0.001]), 4),  # Testunequally spaced monotonically decreasing numpy array
                          (np.asarray([0.94, 0.73, 0.8, 0.21]), 4)])  # Test numpy array that increases and decreases
def test_great_arc_coordinates(points_requested, points_expected):

    m = sunpy.map.Map(sunpy.map.Map(sunpy.data.test.get_test_filepath('aia_171_level1.fits')))
    coordinate_frame = m.coordinate_frame
    a = SkyCoord(600*u.arcsec, -600*u.arcsec, frame=coordinate_frame)
    b = SkyCoord(-100*u.arcsec, 800*u.arcsec, frame=coordinate_frame)
    gc = GreatArc(a, b, points=points_requested)
    coordinates = gc.coordinates()

    # Test the calculation of the SkyCoords
    # Ensure a GreatArc object is returned
    assert isinstance(gc, GreatArc)

    # Test the properties of the GreatArc object
    assert gc.start == a
    assert gc.end == b
    assert gc.distance_unit == u.km
    assert gc.observer == a.observer
    assert gc.center.x == 0 * u.km
    assert gc.center.y == 0 * u.km
    assert gc.center.z == 0 * u.km
    print(gc.start_cartesian)

    np.testing.assert_almost_equal(gc.start_cartesian, np.asarray([432287.18770857, -432289.01663457,  331640.42501877]))
    np.testing.assert_almost_equal(gc.end_cartesian, np.asarray([-72022.85263048,  576185.73277392,  382823.83430374]))
    np.testing.assert_almost_equal(gc.center_cartesian == np.asarray([0, 0, 0]))

    # Great arc properties calculation
    # Vector from center to first point
    np.testing.assert_almost_equal(gc.v1, np.asarray([432287.18770857, -432289.01663457,  331640.42501877]))
    np.testing.assert_almost_equal(gc._r, 695508.00000179256)
    np.testing.assert_almost_equal(gc.v2, np.asarray([-72022.85263048,  576185.73277392,  382823.83430374]))
    np.testing.assert_almost_equal(gc.v3, np.asarray([68458.80553439,  463084.74994859,  514390.2063379]))

    # Inner angle
    assert gc.inner_angle.unit == u.rad
    np.testing.assert_almost_equal(gc.inner_angle.value, 1.8931661720012904)

    # Radius of the sphere
    assert gc.radius.unit == u.km
    np.testing.assert_almost_equal(gc.radius.value, 695508.0000017926)

    # Distance on the sphere between the start point and the end point.
    assert gc.distance.unit == u.km
    np.testing.assert_almost_equal(gc.distance, 1316712.217959667)

    """
    # Correct number of points along the great arc.
    assert len(coordinates) == points_expected

    # Start and end point values are as expected
    np.testing.assert_almost_equal(coordinates[0].Tx.value, a.Tx.value)
    np.testing.assert_almost_equal(coordinates[0].Ty.value, a.Ty.value)
    np.testing.assert_almost_equal(coordinates[-1].Tx.value, b.Tx.value)
    np.testing.assert_almost_equal(coordinates[-1].Ty.value, b.Ty.value)

    # Make sure the output observer is correct
    for i in range(0, len(coordinates)):
        assert coordinates[i].observer == a.observer
    """
