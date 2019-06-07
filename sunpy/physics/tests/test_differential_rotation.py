import os
import pytest

import numpy as np
import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.coordinates import Longitude
from astropy.tests.helper import assert_quantity_allclose
from astropy.time import TimeDelta

from sunpy.coordinates import frames
from sunpy.coordinates.ephemeris import get_earth
from sunpy.map.maputils import map_edges
from sunpy.physics.differential_rotation import (diff_rot, solar_rotate_coordinate,
                                                 differential_rotate,
                                                 _get_new_observer, _rotate_submap_edge,
                                                 _get_extreme_position, _get_bounding_coordinates,
                                                 _warp_sun_coordinates)
import sunpy.data.test
import sunpy.map

# pylint: disable=C0103,R0904,W0201,W0212,W0232,E1103

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


testpath = sunpy.data.test.rootdir

@pytest.fixture
def aia171_test_map():
    return sunpy.map.Map((os.path.join(testpath, 'aia_171_level1.fits')))


@pytest.fixture
def all_off_disk_map(aia171_test_map):
    return aia171_test_map.submap((1, 1)*u.pix, (11, 12)*u.pix)


@pytest.fixture
def all_on_disk_map(aia171_test_map):
    return aia171_test_map.submap((30, 60)*u.pix, (50, 85)*u.pix)


@pytest.fixture
def straddles_limb_map(aia171_test_map):
    return aia171_test_map.submap((64, 80)*u.pix, (120, 127)*u.pix)


@pytest.fixture
def aia171_test_map_with_mask(aia171_test_map):
    shape = aia171_test_map.data.shape
    mask = np.zeros_like(aia171_test_map.data, dtype=bool)
    mask[0:shape[0]//2, 0:shape[1]//2] = True
    return sunpy.map.Map(np.ma.array(aia171_test_map.data, mask=mask), aia171_test_map.meta)


@pytest.fixture
def aia171_test_submap(aia171_test_map):
    bl = SkyCoord(-512 * u.arcsec,  100 * u.arcsec, frame=aia171_test_map.coordinate_frame)
    ur = SkyCoord(-100 * u.arcsec, 400 * u.arcsec, frame=aia171_test_map.coordinate_frame)
    return aia171_test_map.submap(bl, ur)


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
    obs_time = '2010-09-10 12:34:56'
    observer = get_earth(obs_time)
    c = SkyCoord(-570*u.arcsec, 120*u.arcsec, obstime=obs_time, observer=observer, frame=frames.Helioprojective)
    new_time = '2010-09-11 12:34:56'
    new_observer = get_earth(new_time)

    # Test that when both the observer and the time are specified, an error is raised.
    with pytest.raises(ValueError):
        d = solar_rotate_coordinate(c, observer=observer, time=new_time)

    # Test that the code properly filters the observer keyword
    with pytest.raises(ValueError):
        d = solar_rotate_coordinate(c, observer='earth')

    # Test that the code properly filters the time keyword
    with pytest.raises(ValueError):
        d = solar_rotate_coordinate(c, time='noon')

    # Test that the code gives the same output for multiple different inputs
    # that define the same observer location and time.
    for i, definition in enumerate((1 * u.day, TimeDelta(1*u.day), new_time, new_observer)):
        if i in (0, 1, 2):
            d = solar_rotate_coordinate(c, time=definition)
        else:
            d = solar_rotate_coordinate(c, observer=definition)

        # Test that a SkyCoordinate is created
        assert isinstance(d, SkyCoord)

        # Test the coordinate
        np.testing.assert_almost_equal(d.Tx.to(u.arcsec).value, -371.8885208634674, decimal=1)
        np.testing.assert_almost_equal(d.Ty.to(u.arcsec).value, 105.35006656251727, decimal=1)
        np.testing.assert_allclose(d.distance.to(u.km).value, 1.499642e+08, rtol=1e-5)

        # Test that the SkyCoordinate is Helioprojective
        assert isinstance(d.frame, frames.Helioprojective)


def test_differential_rotate(aia171_test_map, all_off_disk_map, all_on_disk_map, straddles_limb_map):

    # Test a map that is entirely off the disk of the Sun
    # Should report an error
    with pytest.raises(ValueError):
        dmap = differential_rotate(all_off_disk_map)

    # Test a full disk map
    new_observer = get_earth(aia171_test_map.date + 6*u.hr)
    dmap = differential_rotate(aia171_test_map, observer=new_observer)
    assert dmap.data.shape == aia171_test_map.data.shape

    # Test a map that is entirely on disk - triggers sub full disk branches
    # Rotated map should have a smaller extent in the x - direction
    new_observer = get_earth(all_on_disk_map.date - 48*u.hr)
    dmap = differential_rotate(all_on_disk_map, observer=new_observer)
    assert dmap.data.shape[1] < all_on_disk_map.data.shape[1]

    # This rotated map should have a larger extent in the x direction
    new_observer = get_earth(all_on_disk_map.date + 48*u.hr)
    dmap = differential_rotate(all_on_disk_map, observer=new_observer)
    assert dmap.data.shape[1] > all_on_disk_map.data.shape[1]

    # Test a map that straddles the limb - triggers sub full disk branches
    # Rotated map should have a smaller extent in the x - direction
    new_observer = get_earth(straddles_limb_map.date + 48*u.hr)
    dmap = differential_rotate(straddles_limb_map, observer=new_observer)
    assert dmap.data.shape[1] < straddles_limb_map.data.shape[1]

    # The output map should have the positional properties of the observer
    assert dmap.date == new_observer.obstime
    assert dmap.heliographic_latitude == new_observer.lat
    assert dmap.heliographic_longitude == new_observer.lon


# Tests of the helper functions
def test_get_new_observer(aia171_test_map):
    initial_obstime = aia171_test_map.date
    rotation_interval = 2 * u.day
    new_time = initial_obstime + rotation_interval
    time_delta = new_time - initial_obstime
    observer = get_earth(initial_obstime + rotation_interval)

    # The observer time is set along with other definitions of time
    for time in (rotation_interval, new_time, time_delta):
        with pytest.raises(ValueError):
            new_observer = _get_new_observer(initial_obstime, observer, time)

    # Obstime property is present but the value is None
    observer_obstime_is_none = SkyCoord(12*u.deg, 46*u.deg, frame=frames.HeliographicStonyhurst)
    with pytest.raises(ValueError):
        new_observer = _get_new_observer(None, observer_obstime_is_none, None)

    # When the observer is set, it gets passed back out
    new_observer = _get_new_observer(initial_obstime, observer, None)
    assert isinstance(new_observer, SkyCoord)
    np.testing.assert_almost_equal(new_observer.transform_to(frames.HeliographicStonyhurst).lon.to(u.deg).value,
                                   observer.transform_to(frames.HeliographicStonyhurst).lon.to(u.deg).value, decimal=3)
    np.testing.assert_almost_equal(new_observer.transform_to(frames.HeliographicStonyhurst).lat.to(u.deg).value,
                                   observer.transform_to(frames.HeliographicStonyhurst).lat.to(u.deg).value, decimal=3)
    np.testing.assert_almost_equal(new_observer.transform_to(frames.HeliographicStonyhurst).radius.to(u.au).value,
                                   observer.transform_to(frames.HeliographicStonyhurst).radius.to(u.au).value, decimal=3)

    # When the time is set, a coordinate for Earth comes back out
    for time in (rotation_interval, new_time, time_delta):
        new_observer = _get_new_observer(initial_obstime, None, time)
        assert isinstance(new_observer, SkyCoord)

        np.testing.assert_almost_equal(new_observer.transform_to(frames.HeliographicStonyhurst).lon.to(u.deg).value,
                                       observer.transform_to(frames.HeliographicStonyhurst).lon.to(u.deg).value, decimal=3)
        np.testing.assert_almost_equal(new_observer.transform_to(frames.HeliographicStonyhurst).lat.to(u.deg).value,
                                       observer.transform_to(frames.HeliographicStonyhurst).lat.to(u.deg).value, decimal=3)
        np.testing.assert_almost_equal(new_observer.transform_to(frames.HeliographicStonyhurst).radius.to(u.au).value,
                                       observer.transform_to(frames.HeliographicStonyhurst).radius.to(u.au).value, decimal=3)

    # The observer and the time cannot both be None
    with pytest.raises(ValueError):
        new_observer = _get_new_observer(initial_obstime, None, None)


def test_rotate_submap_edge(aia171_test_map, all_off_disk_map, all_on_disk_map, straddles_limb_map):

    observer = get_earth(aia171_test_map.date + 2*u.day)

    # For a map that has all the edges off disk, the function should
    # return just the edges of the map - no solar rotation applied.
    for this_map in (aia171_test_map, all_off_disk_map):
        edges = map_edges(this_map)
        for this_edge in range(0, 4):
            pixels = edges[this_edge]
            res = _rotate_submap_edge(this_map, pixels, observer)
            assert all(res.Tx == (this_map.pixel_to_world(pixels[:, 0], pixels[:, 1])).Tx)
            assert all(res.Ty == (this_map.pixel_to_world(pixels[:, 0], pixels[:, 1])).Ty)

    # For an on disk map, all the edges should change
    edges = map_edges(all_on_disk_map)
    for this_edge in range(0, 4):
        pixels = edges[this_edge]
        res = _rotate_submap_edge(all_on_disk_map, pixels, observer)
        assert all(res.Tx != (all_on_disk_map.pixel_to_world(pixels[:, 0], pixels[:, 1])).Tx)
        assert all(res.Ty != (all_on_disk_map.pixel_to_world(pixels[:, 0], pixels[:, 1])).Ty)

    # For the limb map, two of the edges move and two do not
    edges = map_edges(straddles_limb_map)
    for this_edge in (0, 3):  # Top and right edges do not move
        pixels = edges[this_edge]
        res = _rotate_submap_edge(straddles_limb_map, pixels, observer)
        assert all(res.Tx == (straddles_limb_map.pixel_to_world(pixels[:, 0], pixels[:, 1])).Tx)
        assert all(res.Ty == (straddles_limb_map.pixel_to_world(pixels[:, 0], pixels[:, 1])).Ty)

    for this_edge in (1, 2):  # Bottom and left edges do move
        pixels = edges[this_edge]
        res = _rotate_submap_edge(straddles_limb_map, pixels, observer)
        assert all(res.Tx != (straddles_limb_map.pixel_to_world(pixels[:, 0], pixels[:, 1])).Tx)
        assert all(res.Ty != (straddles_limb_map.pixel_to_world(pixels[:, 0], pixels[:, 1])).Ty)


def test_get_extreme_position():
    coords = SkyCoord([-1, 0, 1, np.nan]*u.arcsec, [-2, 0, 2, -np.nan]*u.arcsec, frame=frames.Helioprojective)

    assert _get_extreme_position(coords, 'Tx', operator=np.nanmin) == -1
    assert _get_extreme_position(coords, 'Ty', operator=np.nanmin) == -2

    assert _get_extreme_position(coords, 'Tx', operator=np.nanmax) == 1
    assert _get_extreme_position(coords, 'Ty', operator=np.nanmax) == 2

    with pytest.raises(ValueError):
        _get_extreme_position(coords, 'lon', operator=np.nanmax)


def test_get_bounding_coordinates():
    coords = SkyCoord([-1, 0, 1] * u.arcsec, [-2, 0, 2] * u.arcsec, frame=frames.Helioprojective,
                      observer=get_earth("1999-09-13 00:00:00"))
    bl, tr = _get_bounding_coordinates(coords)

    assert bl.Tx == -1*u.arcsec
    assert bl.Ty == -2*u.arcsec
    assert bl.observer == coords[0].observer

    assert tr.Tx == 1*u.arcsec
    assert tr.Ty == 2*u.arcsec
    assert tr.observer == coords[0].observer


def test_warp_sun_coordinates(all_on_disk_map):
    # Define an observer
    new_observer = get_earth(all_on_disk_map.date + 6*u.hr)

    dummy_array = np.zeros((500, 2))

    # Call the warp
    xy2 = _warp_sun_coordinates(dummy_array, all_on_disk_map, new_observer)

    # Test the properties of the output
    shape = all_on_disk_map.data.shape
    assert xy2.shape == dummy_array.shape
    assert isinstance(xy2, np.ndarray)

    # Test the values - values are not independently found
    # We are passing in 500 pairs of (0,0) so all the output pixels should be the same
    np.testing.assert_almost_equal(xy2[:, 0], 2.3776536242111312, decimal=2)
    np.testing.assert_almost_equal(xy2[:, 1], 0.23785889303423113, decimal=2)
