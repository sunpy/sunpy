import os
import pytest

import numpy as np
import astropy.units as u
from astropy.coordinates import SkyCoord, BaseCoordinateFrame

from sunpy.map.maputils import (all_pixel_indices_from_map,
                                all_coordinates_from_map, map_edges, solar_angular_radius,
                                contains_full_disk, is_all_off_disk,
                                is_all_on_disk, contains_limb,
                                coordinate_is_on_solar_disk, on_disk_bounding_coordinates)
from sunpy.coordinates import HeliographicStonyhurst
import sunpy.data.test
import sunpy.map

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


def test_map_edges(all_off_disk_map):
    edges = map_edges(all_off_disk_map)
    assert type(edges) is tuple
    assert len(edges[2]) == 11
    assert np.all(edges[2][0] == [0, 0] * u.pix)
    assert np.all(edges[2][10] == [0, 10] * u.pix)

    assert len(edges[3]) == 11
    assert np.all(edges[3][0] == [9, 0] * u.pix)
    assert np.all(edges[3][10] == [9, 10] * u.pix)

    assert len(edges[1]) == 10
    assert np.all(edges[1][0] == [0, 0] * u.pix)
    assert np.all(edges[1][9] == [9, 0] * u.pix)

    assert len(edges[0]) == 10
    assert np.all(edges[0][0] == [0, 10] * u.pix)
    assert np.all(edges[0][9] == [9, 10] * u.pix)


def test_solar_angular_radius(aia171_test_map):
    on_disk = aia171_test_map.center
    sar = solar_angular_radius(on_disk)
    assert isinstance(sar, u.Quantity)
    np.testing.assert_almost_equal(sar.to(u.arcsec).value, 971.80181131, decimal=1)


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


def test_coordinate_is_on_solar_disk(aia171_test_map, all_off_disk_map, all_on_disk_map, straddles_limb_map):
    off_disk = aia171_test_map.bottom_left_coord
    on_disk = aia171_test_map.center

    # Check for individual coordinates
    assert coordinate_is_on_solar_disk(on_disk)
    assert ~coordinate_is_on_solar_disk(off_disk)

    # Raise the error
    with pytest.raises(ValueError):
        coordinate_is_on_solar_disk(on_disk.transform_to(HeliographicStonyhurst))

    # Check for sets of coordinates
    assert np.any(coordinate_is_on_solar_disk(all_coordinates_from_map(aia171_test_map)))
    assert np.any(~coordinate_is_on_solar_disk(all_coordinates_from_map(aia171_test_map)))
    assert np.all(~coordinate_is_on_solar_disk(all_coordinates_from_map(all_off_disk_map)))
    assert np.all(coordinate_is_on_solar_disk(all_coordinates_from_map(all_on_disk_map)))
    assert np.any(coordinate_is_on_solar_disk(all_coordinates_from_map(straddles_limb_map)))
    assert np.any(~coordinate_is_on_solar_disk(all_coordinates_from_map(straddles_limb_map)))


# Testing values are derived from running the code, not from external sources
def test_on_disk_bounding_coordinates(aia171_test_map):
    bl, tr = on_disk_bounding_coordinates(aia171_test_map)
    np.testing.assert_almost_equal(bl.Tx.to(u.arcsec).value, -954.17124289, decimal=1)
    np.testing.assert_almost_equal(bl.Ty.to(u.arcsec).value, -965.93063472, decimal=1)
    np.testing.assert_almost_equal(tr.Tx.to(u.arcsec).value, 964.27061417, decimal=1)
    np.testing.assert_almost_equal(tr.Ty.to(u.arcsec).value, 971.63586861, decimal=1)
