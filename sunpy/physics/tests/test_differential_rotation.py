from __future__ import absolute_import
import pytest
import numpy as np
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.coordinates import Longitude, Latitude, Angle
from astropy.tests.helper import assert_quantity_allclose
from sunpy.coordinates import frames
from sunpy.coordinates.ephemeris import get_earth
from sunpy.physics.differential_rotation import diff_rot, solar_rotate_coordinate

#pylint: disable=C0103,R0904,W0201,W0212,W0232,E1103

# Please note the numbers in these tests are not checked for physical
# accuracy, only that they are the values the function was outputting upon
# implementation.

@pytest.fixture
def seconds_per_day():
    return 24 * 60 * 60.0 * u.s


def test_single(seconds_per_day):
    rot = diff_rot(10 * seconds_per_day, 30 * u.deg)
    assert_quantity_allclose(rot, 136.8216 * u.deg, rtol=1e-3)


def test_array(seconds_per_day):
    rot = diff_rot(10 * seconds_per_day, np.linspace(-70, 70, 2) * u.deg)
    assert_quantity_allclose(rot, Longitude(np.array([110.2725,  110.2725]) * u.deg), rtol=1e-3)


def test_synodic(seconds_per_day):
    rot = diff_rot(10 * seconds_per_day, 30 * u.deg, rot_type='howard', frame_time='synodic')
    assert_quantity_allclose(rot, 126.9656 * u.deg, rtol=1e-3)


def test_sidereal(seconds_per_day):
    rot = diff_rot(10 * seconds_per_day, 30 * u.deg, rot_type='howard', frame_time='sidereal')
    assert_quantity_allclose(rot, 136.8216 * u.deg, rtol=1e-3)


def test_howard(seconds_per_day):
    rot = diff_rot(10 * seconds_per_day, 30 * u.deg, rot_type='howard')
    assert_quantity_allclose(rot, 136.8216 * u.deg, rtol=1e-3)


def test_allen(seconds_per_day):
    rot = diff_rot(10 * seconds_per_day, 30 * u.deg, rot_type='allen')
    assert_quantity_allclose(rot, 136.9 * u.deg, rtol=1e-3)


def test_snodgrass(seconds_per_day):
    rot = diff_rot(10 * seconds_per_day, 30 * u.deg, rot_type='snodgrass')
    assert_quantity_allclose(rot, 135.4232 * u.deg, rtol=1e-3)


def test_fail(seconds_per_day):
    with pytest.raises(ValueError):
        rot = diff_rot(10 * seconds_per_day, 30 * u.deg, rot_type='garbage')


def test_solar_rotate_coordinate:
    # Testing along the Sun-Earth line, observer is on the Earth
    obstime = '2010-09-10 12:34:56'
    c = SkyCoord(-570*u.arcsec, 120*u.arcsec, obstime=obstime, observer=get_earth(obstime), frame=frames.Helioprojective)
    d = solar_rotate_coordinate(c, '2010-09-10 13:34:56')

    # Test that a SkyCoordinate is created
    assert isinstance(d, SkyCoord)


    <SkyCoord (Helioprojective: obstime=2010-09-10 13:34:56, rsun=695508.0 km, observer=<HeliographicStonyhurst Coordinate (obstime=2010-09-10 13:34:56): (lon, lat, radius) in (deg, deg, AU)
    ( 0.,  7.24822784,  1.00695436)>): (Tx, Ty, distance) in (arcsec, arcsec, km)
    (-562.37689548,  119.26840368,   1.50083152e+08)>
                   '2012-06-15', '2012-06-15 16:05:23')
    np.testing.assert_almost_equal(x.to(u.arcsec).value, 574.2, decimal=1)
    np.testing.assert_almost_equal(y.to(u.arcsec).value, -108.4, decimal=1)
    # Test that astropy Angles are returned and that they have the expected
    # units

