import os
import pytest

import numpy as np
from astropy import units as u
from astropy.coordinates import SkyCoord, BaseCoordinateFrame
from astropy.coordinates import Longitude
from astropy.tests.helper import assert_quantity_allclose
from astropy.time import TimeDelta

from sunpy.coordinates import frames
from sunpy.coordinates.ephemeris import get_earth
from sunpy.physics.differential_rotation import diff_rot, solar_rotate_coordinate,\
    diffrot_map, all_pixel_indices_from_map, all_coordinates_from_map,\
    find_pixel_radii, map_edges, contains_full_disk, is_all_off_disk,\
    is_all_on_disk, contains_limb, coordinate_is_on_disk, on_disk_bounding_coordinates,\
    _get_new_observer, _rotate_submap_edge, _get_extreme_position, _get_bounding_coordinates,\
    _warp_sun_coordinates
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
    return aia171_test_map.submap((50, 60)*u.pix, (70, 85)*u.pix)


@pytest.fixture
def straddles_limb_map(aia171_test_map):
    return aia171_test_map.submap((0, 0)*u.pix, (50, 40)*u.pix)


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
    new_time = '2010-09-10 13:34:56'
    new_observer = get_earth(new_time)

    # Test that when both the observer and the time are specified, an error is raised.
    with pytest.raises(ValueError):
        d = solar_rotate_coordinate(c, observer=observer, time=new_time)

    # Test that the code properly filters the observer keyword
    with pytest.raises(ValueError):
        d = solar_rotate_coordinate(c, observer='earth')
    with pytest.raises(TypeError):
        no_obstime = SkyCoord(-570*u.arcsec, 120*u.arcsec, frame=frames.Helioprojective)
        d = solar_rotate_coordinate(c, observer=no_obstime)

    # Test that the code properly filters the time keyword
    with pytest.raises(ValueError):
        d = solar_rotate_coordinate(c, time='noon')

    # Test that the code gives the same output for multiple different inputs
    # that define the same observer location and time.
    for i, definition in enumerate((1 * u.hour, TimeDelta(1*u.hour), new_time, new_observer)):
        if i in (0, 1, 2):
            d = solar_rotate_coordinate(c, time=definition)
        else:
            d = solar_rotate_coordinate(c, observer=definition)

        # Test that a SkyCoordinate is created
        assert isinstance(d, SkyCoord)

        # Test the coordinate
        np.testing.assert_almost_equal(d.Tx.to(u.arcsec).value, -562.89877818, decimal=1)
        np.testing.assert_almost_equal(d.Ty.to(u.arcsec).value, 119.3152842, decimal=1)
        np.testing.assert_almost_equal(d.distance.to(u.km).value, 1.500850782e+08, decimal=1)

        # Test that the SkyCoordinate is Helioprojective
        assert isinstance(d.frame, frames.Helioprojective)


def test_diffrot_map(aia171_test_map, all_off_disk_map, all_on_disk_map, straddles_limb_map):

    # Test that the function correctly identifies a map that is entirely off
    # the disk of the Sun and reports an error correctly
    with pytest.raises(ValueError):
        dmap = diffrot_map(all_off_disk_map)

    with pytest.raises(ValueError):
        dmap =

    # TODO - a lot more tests here, maybe use image hash comparison?

    aia171_test_map_diff_rot =


@pytest.fixture
def sub_smap(aia171_test_map):
    return aia171_test_map.submap((0, 0)*u.pix, (50, 60)*u.pix)


def test_all_pixel_indices_from_map(sub_smap):
    pixel_indices = all_pixel_indices_from_map(sub_smap)
    shape = sub_smap.data.shape
    ny = shape[0]
    nx = shape[1]
    assert np.all(pixel_indices.shape == (2, ny, nx))
    assert np.all(pixel_indices.unit == u.pix)
    assert np.all(pixel_indices[:, 0, 0] == [0., 0.] * u.pix)
    assert np.all(pixel_indices[:, 0, nx-1] == [nx-1, 0.] * u.pix)
    assert np.all(pixel_indices[:, ny-1, 0] == [0., ny-1] * u.pix)
    assert np.all(pixel_indices[:, ny-1, nx-1] == [nx-1, ny-1] * u.pix)


def test_all_coordinates_from_map(sub_smap):
    coordinates = all_coordinates_from_map(sub_smap)
    shape = sub_smap.data.shape
    assert coordinates.shape == (shape[0], shape[1])
    assert isinstance(coordinates, SkyCoord)
    assert isinstance(coordinates.frame, BaseCoordinateFrame)
    assert coordinates.frame.name == sub_smap.coordinate_frame.name


def test_find_pixel_radii(aia171_test_map):
    # The known maximum radius
    known_maximum_pixel_radius = 1.77805631 * u.R_sun

    # Calculate the pixel radii
    pixel_radii = find_pixel_radii(aia171_test_map)

    # The shape of the pixel radii is the same as the input map
    assert pixel_radii.shape[0] == int(aia171_test_map.dimensions[0].value)
    assert pixel_radii.shape[1] == int(aia171_test_map.dimensions[1].value)

    # Make sure the unit is solar radii
    assert pixel_radii.unit == u.R_sun

    # Make sure the maximum
    assert np.allclose(np.max(pixel_radii).value, known_maximum_pixel_radius.value)

    # Test that the new scale is used
    pixel_radii = find_pixel_radii(aia171_test_map, scale=2*aia171_test_map.rsun_obs)
    assert np.allclose(np.max(pixel_radii).value, known_maximum_pixel_radius.value / 2)


def test_map_edges(all_off_disk_map):
    edges = map_edges(all_off_disk_map)
    assert type(edges) is dict
    keys = edges.keys()
    assert 'lhs' in keys
    assert 'rhs' in keys
    assert 'top' in keys
    assert 'bottom' in keys
    assert len(edges['lhs']) == 11
    assert np.all(edges['lhs'][0] == [0, 0] * u.pix)
    assert np.all(edges['lhs'][10] == [10, 0] * u.pix)

    assert len(edges['rhs']) == 11
    assert np.all(edges['rhs'][0] == [0, 9] * u.pix)
    assert np.all(edges['rhs'][10] == [10, 9] * u.pix)

    assert len(edges['top']) == 10
    assert np.all(edges['top'][0] == [0, 0] * u.pix)
    assert np.all(edges['top'][9] == [0, 9] * u.pix)

    assert len(edges['bottom']) == 10
    assert np.all(edges['bottom'][0] == [10, 0] * u.pix)
    assert np.all(edges['bottom'][9] == [10, 9] * u.pix)


def test_contains_full_disk(aia171_test_map, all_off_disk_map, all_on_disk_map, straddles_limb_map):
    assert contains_full_disk(aia171_test_map)
    assert ~contains_full_disk(all_off_disk_map)
    assert ~contains_full_disk(all_on_disk_map)
    assert ~contains_full_disk(straddles_limb_map)


def test_is_all_off_disk(aia171_test_map, all_off_disk_map, all_on_disk_map, straddles_limb_map):
    assert ~is_all_off_disk(aia171_test_map)
    assert is_all_off_disk(all_off_disk_map)
    assert ~is_all_off_disk(all_on_disk_map)
    assert ~is_all_off_disk(straddles_limb_map)


def test_is_all_on_disk(aia171_test_map, all_off_disk_map, all_on_disk_map, straddles_limb_map):
    assert ~is_all_on_disk(aia171_test_map)
    assert ~is_all_on_disk(all_off_disk_map)
    assert is_all_on_disk(all_on_disk_map)
    assert ~is_all_on_disk(straddles_limb_map)


def test_contains_limb(aia171_test_map, all_off_disk_map, all_on_disk_map, straddles_limb_map):
    assert contains_limb(aia171_test_map)
    assert ~contains_limb(all_off_disk_map)
    assert ~contains_limb(all_on_disk_map)
    assert contains_limb(straddles_limb_map)


def test_coordinate_is_on_disk(aia171_test_map, all_off_disk_map, all_on_disk_map, straddles_limb_map):
    scale = aia171_test_map.rsun_obs
    off_disk = aia171_test_map.bottom_left_coord
    on_disk = aia171_test_map.center

    # Check for individual coordinates
    assert coordinate_is_on_disk(on_disk, scale)
    assert ~coordinate_is_on_disk(off_disk, scale)

    # Check for sets of coordinates
    assert np.any(coordinate_is_on_disk(all_coordinates_from_map(aia171_test_map), scale))
    assert np.any(~coordinate_is_on_disk(all_coordinates_from_map(aia171_test_map), scale))
    assert np.all(~coordinate_is_on_disk(all_coordinates_from_map(all_off_disk_map), scale))
    assert np.all(coordinate_is_on_disk(all_coordinates_from_map(all_on_disk_map), scale))
    assert np.any(coordinate_is_on_disk(all_coordinates_from_map(straddles_limb_map), scale))
    assert np.any(~coordinate_is_on_disk(all_coordinates_from_map(straddles_limb_map), scale))


# Testing values are derived from running the code, not from external sources
def test_on_disk_bounding_coordinates(aia171_test_map):
    bl, tr = on_disk_bounding_coordinates(aia171_test_map)
    np.testing.assert_almost_equal(bl.Tx.to(u.arcsec).value, -954.17124289, decimal=1)
    np.testing.assert_almost_equal(bl.Ty.to(u.arcsec).value, -965.93063472, decimal=1)
    np.testing.assert_almost_equal(tr.Tx.to(u.arcsec).value, 964.27061417, decimal=1)
    np.testing.assert_almost_equal(tr.Ty.to(u.arcsec).value, 971.63586861, decimal=1)


# Tests of the helper functions
def test_get_new_observer(aia171_test_map):
    initial_obstime = aia171_test_map.date
    rotation_interval = 2 * u.day
    new_time = initial_obstime + rotation_interval
    time_delta = new_time - initial_obstime
    observer = get_earth(initial_obstime + rotation_interval)

    # The observer and the time cannot both be not None
    for time in (rotation_interval, new_time, time_delta):
        with pytest.raises(ValueError):
            new_observer = _get_new_observer(initial_obstime, observer, time)

    # When the observer is set, it gets passed back out
    new_observer = _get_new_observer(initial_obstime, observer, None)
    assert isinstance(new_observer, SkyCoord)
    assert new_observer.transform_to(frames.HeliographicStonyhurst).lon == observer.transform_to(frames.HeliographicStonyhurst).lon
    assert new_observer.transform_to(frames.HeliographicStonyhurst).lat == observer.transform_to(frames.HeliographicStonyhurst).lat
    assert new_observer.transform_to(frames.HeliographicStonyhurst).radius == observer.transform_to(frames.HeliographicStonyhurst).radius

    # When the time is set, a coordinate for Earth comes back out
    for time in (rotation_interval, new_time, time_delta):
        new_observer = _get_new_observer(initial_obstime, None, time)
        assert isinstance(new_observer, SkyCoord)
        assert new_observer.transform_to(frames.HeliographicStonyhurst).lon == observer.transform_to(
            frames.HeliographicStonyhurst).lon
        assert new_observer.transform_to(frames.HeliographicStonyhurst).lat == observer.transform_to(
            frames.HeliographicStonyhurst).lat
        assert new_observer.transform_to(frames.HeliographicStonyhurst).radius == observer.transform_to(
            frames.HeliographicStonyhurst).radius

    # The observer and the time cannot both be None
    with pytest.raises(ValueError):
        new_observer = _get_new_observer(initial_obstime, None, None)


def test_rotate_submap_edge(aia171_test_map, all_off_disk_map, all_on_disk_map, straddles_limb_map):

    observer = get_earth(aia171_test_map.date + 2*u.day)

    # For a map that has all the edges off disk, the function should
    # return just the edges of the map - no solar rotation applied.
    for this_map in (aia171_test_map, all_off_disk_map):
        edges = map_edges(this_map)
        for this_edge in edges.keys():
            pixels = edges[this_edge]
            res = _rotate_submap_edge(this_map, pixels, observer)
            assert all(res.Tx == (this_map.pixel_to_world(pixels[:, 1], pixels[:, 0])).Tx)
            assert all(res.Ty == (this_map.pixel_to_world(pixels[:, 1], pixels[:, 0])).Ty)

    # For an on disk map, all the edges should change
    edges = map_edges(all_on_disk_map)
    for this_edge in edges.keys():
        pixels = edges[this_edge]
        res = _rotate_submap_edge(all_on_disk_map, pixels, observer)
        assert all(res.Tx != (all_on_disk_map.pixel_to_world(pixels[:, 1], pixels[:, 0])).Tx)
        assert all(res.Ty != (all_on_disk_map.pixel_to_world(pixels[:, 1], pixels[:, 0])).Ty)

    # For the limb map, two of the edges move and two do not
    edges = map_edges(straddles_limb_map)
    for this_edge in ('top', 'rhs'):
        pixels = edges[this_edge]
        res = _rotate_submap_edge(all_on_disk_map, pixels, observer)
        assert all(res.Tx != (all_on_disk_map.pixel_to_world(pixels[:, 1], pixels[:, 0])).Tx)
        assert all(res.Ty != (all_on_disk_map.pixel_to_world(pixels[:, 1], pixels[:, 0])).Ty)

    for this_edge in ('bottom', 'lhs'):
        pixels = edges[this_edge]
        res = _rotate_submap_edge(all_on_disk_map, pixels, observer)
        assert all(res.Tx != (all_on_disk_map.pixel_to_world(pixels[:, 1], pixels[:, 0])).Tx)
        assert all(res.Ty != (all_on_disk_map.pixel_to_world(pixels[:, 1], pixels[:, 0])).Ty)


def test_get_extreme_position():
    coords = SkyCoord([-1, 0, 1, np.nan]*u.arcsec, [-2, 0, 2, -np.nan]*u.arcsec, frame=frames.Helioprojective)

    assert _get_extreme_position(coords, 'Tx', operator=np.nanmin) == -1*u.arcsec
    assert _get_extreme_position(coords, 'Ty', operator=np.nanmin) == -2*u.arcsec

    assert _get_extreme_position(coords, 'Tx', operator=np.nanmax) == 1*u.arcsec
    assert _get_extreme_position(coords, 'Ty', operator=np.nanmax) == 2*u.arcsec

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
    assert tr.Tt == 2*u.arcsec
    assert tr.observer == coords[0].observer


def test_warp_sun_coordinates(all_on_disk_map):
    new_observer = get_earth(all_on_disk_map.date + 1*u.hr)

    # This array is not used in the function but is part of the signature
    dummy_array = np.zeroes(10)

    # Call the warp
    xy2 = _warp_sun_coordinates(dummy_array, all_on_disk_map, new_observer)

    # Test the properties of the output
    assert xy2.shape == 0  # fill this in with the correct shape
    assert isinstance(xy2, np.masked_array)  # fill this in with the correct type
