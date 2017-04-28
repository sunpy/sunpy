from __future__ import absolute_import
import pytest
import numpy as np
from astropy import units as u
from astropy.coordinates import Longitude, Latitude, Angle
from sunpy.physics.differential_rotation import diff_rot, _sun_pos, _calc_P_B0_SD, rot_hpc
from sunpy.tests.helpers import assert_quantity_allclose
#pylint: disable=C0103,R0904,W0201,W0212,W0232,E1103

# Please note the numbers in these tests are not checked for physical
# accuracy, only that they are the values the function was outputting upon
# implementation.

@pytest.fixture
def seconds_per_day():
    return 24 * 60 * 60.0 * u.s

def test_single(seconds_per_day):
    rot = diff_rot(10 * seconds_per_day, 30 * u.deg)
    assert rot == 136.8216 * u.deg


def test_array(seconds_per_day):
    rot = diff_rot(10 * seconds_per_day, np.linspace(-70, 70, 2) * u.deg)
    assert_quantity_allclose(rot, Longitude(np.array([110.2725,  110.2725]) * u.deg))


def test_synodic(seconds_per_day):
    rot = diff_rot(10 * seconds_per_day, 30 * u.deg, rot_type='howard', frame_time='synodic')
    assert rot == 126.9656 * u.deg


def test_sidereal(seconds_per_day):
    rot = diff_rot(10 * seconds_per_day, 30 * u.deg, rot_type='howard', frame_time='sidereal')
    assert rot == 136.8216 * u.deg


def test_howard(seconds_per_day):
    rot = diff_rot(10 * seconds_per_day, 30 * u.deg, rot_type='howard')
    assert rot == 136.8216 * u.deg


def test_allen(seconds_per_day):
    rot = diff_rot(10 * seconds_per_day, 30 * u.deg, rot_type='allen')
    assert rot == 136.9 * u.deg


def test_snodgrass(seconds_per_day):
    rot = diff_rot(10 * seconds_per_day, 30 * u.deg, rot_type='snodgrass')
    assert rot == 135.4232 * u.deg


def test_fail(seconds_per_day):
    with pytest.raises(ValueError):
        rot = diff_rot(10 * seconds_per_day, 30 * u.deg, rot_type='garbage')


def test_sunpos():
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


def test_calc_P_B0_SD():
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


def test_rot_hpc():
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
