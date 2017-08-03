from __future__ import absolute_import
import pytest
import numpy as np
from astropy import units as u
from astropy.time import Time
from astropy.coordinates import SkyCoord
from astropy.coordinates import Longitude, Latitude, Angle
from astropy.tests.helper import assert_quantity_allclose
from sunpy.coordinates import frames
from sunpy.coordinates.ephemeris import get_earth
from sunpy.physics.differential_rotation import diff_rot, solar_rotate_coordinate
from sunpy.time import parse_time

#pylint: disable=C0103,R0904,W0201,W0212,W0232,E1103

# Please note the numbers in these tests are not checked for physical
# accuracy, only that they are the values the function was outputting upon
# implementation.  This is not a significant issue for the diff_rot function
# since it is relatively simple and the values it produces can be easily
# compared to other implementations of the same simple function.  The same
# cannot be said for the solar_rotate_coordinate function.  This functionality
# relies accurate knowledge of the solar ephemeris in particular.
# There is no reference implementation of the solar_rotate_coordinate function
# of demonstrated trustworthiness at time of writing in any language.  There
# are no known independent values or tests that can be used to test the
# veracity of the solar_rotate_coordinate function.  This being the case, the
# solar_rotate_coordinate function is tested against values that it generated.
# Therefore these tests test for consistency, not accuracy.  Note that when the
# 0.8.0 branch was released, the solar ephemeris calculation was handed off to
# the relevant Astropy code.  The solar_rotate_coordinate tests were changed
# for self-consistency.  Note that the change in position comparing the results
# of pre- and 0.8.0 sunpy solar coordinate rotation functionality (rot_hpc
# and solar_rotate_coordinate respectively) was on the order of 0.5 arcseconds.
# At time of writing, the difference between the rotation
# calculated using the pre-0.8.0 rot_hpc function and the SSWIDL equivalent
# rot_xy.pro for the tests given in pre-0.8.0 were on the order of hundredths
# of an arcsecond.  I suspect that the reason for the small differences is
# because the sunpy's ephemeris and coordinate transformation infrastructure
# was largely based on that in SSWIDL.


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


def test_solar_rotate_coordinate():
    # Testing along the Sun-Earth line, observer is on the Earth
    obstime = '2010-09-10 12:34:56'
    newtime = '2010-09-10 13:34:56'
    c = SkyCoord(-570*u.arcsec, 120*u.arcsec, obstime=obstime, observer=get_earth(obstime), frame=frames.Helioprojective)
    d = solar_rotate_coordinate(c, newtime)

    # Test that a SkyCoordinate is created
    assert isinstance(d, SkyCoord)

    # Test the coordinate
    np.testing.assert_almost_equal(d.Tx.to(u.arcsec).value, -562.3768, decimal=1)
    np.testing.assert_almost_equal(d.Ty.to(u.arcsec).value, 119.2684, decimal=1)
    np.testing.assert_almost_equal(d.distance.to(u.km).value, 150083151.97246578, decimal=1)

    # Test that the SkyCoordinate is Helioprojective
    assert isinstance(d.frame, frames.Helioprojective)

    # Test the observer
    assert d.observer.obstime == Time(parse_time(newtime), scale='utc')
    np.testing.assert_almost_equal(d.observer.lon.to(u.deg).value, 0.0, decimal=5)
    np.testing.assert_almost_equal(d.observer.lat.to(u.deg).value, 7.248, decimal=3)
    np.testing.assert_almost_equal(d.observer.radius.to(u.AU).value, 1.006954, decimal=6)
    assert isinstance(d.observer, frames.HeliographicStonyhurst)
