from __future__ import absolute_import
<<<<<<< HEAD
import pytest
import numpy as np
from astropy import units as u
from astropy.units import Quantity
from astropy.coordinates import Longitude, Latitude, Angle
from sunpy.physics.transforms.differential_rotation import diff_rot, _sun_pos, _calc_P_B0_SD, rot_hpc
#pylint: disable=C0103,R0904,W0201,W0212,W0232,E1103


class TestDiffRot():
    """ Please note the numbers in these tests are not checked for physical
    accuracy, only that they are the values the function was outputting upon
    implementation. """

    def test_single(self):
        rot = diff_rot(10, 30 * u.deg)
        np.testing.assert_almost_equal(rot.to(u.deg).value, 136.8216, decimal=4)
        isinstance(rot, Longitude)
        rot.unit == u.deg

    def test_array(self):
        rot = diff_rot(10, np.linspace(-70, 70, 2) * u.deg)
        np.testing.assert_almost_equal(rot.to(u.deg).value,
                                       np.array([110.2725, 110.2725]),
                                       decimal=4)

    def test_synodic(self):
        rot = diff_rot(10, 30 * u.deg, rot_type='howard', frame_time='synodic')
        np.testing.assert_almost_equal(rot.to(u.deg).value, 126.9656, decimal=4)

    def test_sidereal(self):
        rot = diff_rot(10, 30 * u.deg, rot_type='howard', frame_time='sidereal')
        np.testing.assert_almost_equal(rot.to(u.deg).value, 136.8216, decimal=4)

    def test_howard(self):
        rot = diff_rot(10, 30 * u.deg, rot_type='howard')
        np.testing.assert_almost_equal(rot.to(u.deg).value, 136.8216, decimal=4)

    def test_allen(self):
        rot = diff_rot(10, 30 * u.deg, rot_type='allen')
        np.testing.assert_almost_equal(rot.to(u.deg).value, 136.9, decimal=4)

    def test_snodgrass(self):
        rot = diff_rot(10, 30 * u.deg, rot_type='snodgrass')
        np.testing.assert_almost_equal(rot.to(u.deg).value, 135.4232, decimal=4)

    def test_fail(self):
        pytest.raises(ValueError, diff_rot, 10, 30 * u.deg, rot_type='garbage')

    def test_sunpos(self):
        result = _sun_pos('2013-05-14')
        assertion = {'obliq': (23.4358, Angle, u.deg),
                     'app_long': (53.3683, Longitude, u.deg),
                     'dec': (18.6125, Latitude, u.deg),
                     'ra': (50.9796, Longitude, u.deg),
                     'longitude': (53.3705, Longitude, u.deg)}
        for k in assertion:
            np.testing.assert_almost_equal(result[k].to(u.deg).value, assertion[k][0], decimal=4)
            isinstance(result[k], assertion[k][1])
            result[k].unit == assertion[k][2]

    def test_calc_P_B0_SD(self):
        result = _calc_P_B0_SD('2012-12-14')
        assertion = {'p': (10.4868, Angle, u.deg),
                     'b0': (-0.8127, Angle, u.deg),
                     'l0': (0.0000, Angle, u.deg),
                     'sd': (16.2364 / 60.0, Angle, u.arcmin)}
        for k in assertion:
            np.testing.assert_almost_equal(result[k].to(u.degree).value,
                                           assertion[k][0], decimal=4)
            # Test that the correct astropy Quantity objects are returned and
            # that they have the expected units.
            isinstance(result[k], assertion[k][1])
            result[k].unit == assertion[k][2]

    def test_rot_hpc(self):
        # testing along the Sun-Earth line, observer is on the Earth
        x, y = rot_hpc(451.4 * u.arcsec, -108.9 * u.arcsec,
                       '2012-06-15', '2012-06-15 16:05:23')
        np.testing.assert_almost_equal(x.to(u.arcsec).value, 574.2, decimal=1)
        np.testing.assert_almost_equal(y.to(u.arcsec).value, -108.4, decimal=1)
        # Test that astropy Angles are returned and that they have the expected
        # units
        isinstance(x, Angle)
        x.unit == u.arcsec
        isinstance(y, Angle)
        y.unit == u.arcsec
=======

import pytest

import numpy as np
from astropy import units as u
from astropy.coordinates import Longitude

from sunpy.physics.transforms.differential_rotation import diff_rot
from sunpy.tests.helpers import assert_quantity_allclose

# Please note the numbers in these tests are not checked for physical
# accuracy, only that they are the values the function was outputting upon
# implementation.

def test_single():
    rot = diff_rot(10, 30 * u.deg)
    assert rot == 136.8216 * u.deg

def test_array():
    rot = diff_rot(10, np.linspace(-70, 70, 2) * u.deg)
    assert_quantity_allclose(rot, Longitude(np.array([110.2725,  110.2725]) * u.deg))

def test_synodic():
    rot = diff_rot(10, 30 * u.deg, rot_type='howard', frame_time='synodic')
    assert rot == 126.9656 * u.deg

def test_sidereal():
    rot = diff_rot(10, 30 * u.deg, rot_type='howard', frame_time='sidereal')
    assert rot == 136.8216 * u.deg

def test_howard():
    rot = diff_rot(10, 30 * u.deg, rot_type='howard')
    assert rot == 136.8216 * u.deg

def test_allen():
    rot = diff_rot(10, 30 * u.deg, rot_type='allen')
    assert rot == 136.9 * u.deg

def test_snodgrass():
    rot = diff_rot(10, 30 * u.deg, rot_type='snodgrass')
    assert rot == 135.4232 * u.deg

def test_fail():
    with pytest.raises(ValueError):
        rot = diff_rot(10, 30 * u.deg, rot_type='garbage')
>>>>>>> 931049fa586c026fbba9a23c34d6006d2a41a5ab
