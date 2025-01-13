import re

import numpy as np
import pytest

import astropy.units as u
from astropy.coordinates import BaseCoordinateFrame, SkyCoord
from astropy.tests.helper import assert_quantity_allclose

import sunpy.map
from sunpy.coordinates import HeliographicStonyhurst
from sunpy.coordinates.frames import HeliographicCarrington
from sunpy.coordinates.utils import GreatArc
from sunpy.map.maputils import (
    _verify_coordinate_helioprojective,
    all_coordinates_from_map,
    all_corner_coords_from_map,
    all_pixel_indices_from_map,
    contains_coordinate,
    contains_full_disk,
    contains_limb,
    contains_solar_center,
    coordinate_is_on_solar_disk,
    is_all_off_disk,
    is_all_on_disk,
    map_edges,
    on_disk_bounding_coordinates,
    pixelate_coord_path,
    sample_at_coords,
    solar_angular_radius,
)


@pytest.fixture
def all_off_disk_map(aia171_test_map):
    return aia171_test_map.submap((1, 1)*u.pix, top_right=(11, 12)*u.pix)


@pytest.fixture
def all_on_disk_map(aia171_test_map):
    return aia171_test_map.submap((30, 60)*u.pix, top_right=(50, 85)*u.pix)


@pytest.fixture
def straddles_limb_map(aia171_test_map):
    return aia171_test_map.submap((64, 80)*u.pix, top_right=(120, 127)*u.pix)


@pytest.fixture
def sub_smap(aia171_test_map):
    return aia171_test_map.submap((0, 0)*u.pix, top_right=(50, 60)*u.pix)


@pytest.fixture
def non_helioprojective_map():
    data = np.arange(0, 100).reshape(10, 10)
    coord = SkyCoord(0*u.arcsec, 0*u.arcsec, obstime='2013-10-28 08:24',
                     observer='earth', frame=HeliographicCarrington)
    header = sunpy.map.header_helper.make_fitswcs_header(data, coord,
                                                         reference_pixel=[0, 0]*u.pixel,
                                                         scale=[2, 2]*u.arcsec/u.pixel,
                                                         telescope='Fake Telescope', instrument='UV detector',
                                                         wavelength=1000*u.angstrom)
    return sunpy.map.Map(data, header)


@pytest.fixture
def non_helioprojective_skycoord():
    return SkyCoord(0 * u.rad, 0 * u.rad, frame="icrs")


@pytest.fixture
def aia_test_arc(aia171_test_map):
    start = SkyCoord(735 * u.arcsec, -471 * u.arcsec, frame=aia171_test_map.coordinate_frame)
    end = SkyCoord(-100 * u.arcsec, 800 * u.arcsec, frame=aia171_test_map.coordinate_frame)
    return GreatArc(start, end)


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

    xpix, ypix = sub_smap.world_to_pixel(coordinates[0, 0])
    assert_quantity_allclose(xpix, 0*u.pix, atol=1e-7*u.pix)
    assert_quantity_allclose(ypix, 0*u.pix, atol=1e-7*u.pix)

    xpix, ypix = sub_smap.world_to_pixel(coordinates[-1, -1])
    assert_quantity_allclose(xpix, sub_smap.dimensions[0] - 1*u.pix)
    assert_quantity_allclose(ypix, sub_smap.dimensions[1] - 1*u.pix)


def test_all_corner_coordinates_from_map(sub_smap):
    coordinates = all_corner_coords_from_map(sub_smap)
    shape = sub_smap.data.shape
    assert coordinates.shape == (shape[0] + 1, shape[1] + 1)
    assert isinstance(coordinates, SkyCoord)
    assert isinstance(coordinates.frame, BaseCoordinateFrame)
    assert coordinates.frame.name == sub_smap.coordinate_frame.name

    xpix, ypix = sub_smap.world_to_pixel(coordinates[0, 0])
    assert_quantity_allclose(xpix, -0.5*u.pix)
    assert_quantity_allclose(ypix, -0.5*u.pix)

    xpix, ypix = sub_smap.world_to_pixel(coordinates[-1, -1])
    assert_quantity_allclose(xpix, sub_smap.dimensions[0] - 0.5*u.pix)
    assert_quantity_allclose(ypix, sub_smap.dimensions[1] - 0.5*u.pix)


def test_map_edges(all_off_disk_map):
    edges = map_edges(all_off_disk_map)
    assert isinstance(edges, tuple)
    assert len(edges[2]) == 12
    assert np.all(edges[2][0] == [0, 0] * u.pix)
    assert np.all(edges[2][11] == [0, 11] * u.pix)

    assert len(edges[3]) == 12
    assert np.all(edges[3][0] == [10, 0] * u.pix)
    assert np.all(edges[3][11] == [10, 11] * u.pix)

    assert len(edges[1]) == 11
    assert np.all(edges[1][0] == [0, 0] * u.pix)
    assert np.all(edges[1][10] == [10, 0] * u.pix)

    assert len(edges[0]) == 11
    assert np.all(edges[0][0] == [0, 11] * u.pix)
    assert np.all(edges[0][10] == [10, 11] * u.pix)


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
    with pytest.raises(ValueError, match=re.escape("The input coordinate(s) is of type HeliographicStonyhurst, but must be in the Helioprojective frame.")):
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


def test_data_at_coordinates(aia171_test_map, aia_test_arc):
    data = sample_at_coords(aia171_test_map, aia_test_arc.coordinates())
    pixels = np.asarray(np.rint(
        aia171_test_map.world_to_pixel(aia_test_arc.coordinates())), dtype=int)
    x = pixels[0, :]
    y = pixels[1, :]
    intensity_along_arc = aia171_test_map.data[y, x] * aia171_test_map.unit
    assert_quantity_allclose(data[0], intensity_along_arc[0])
    assert_quantity_allclose(data[-1], intensity_along_arc[-1])


def test_sample_out_of_bounds(aia171_test_map):
    point = aia171_test_map.pixel_to_world([-1, 1]*u.pix, [-1, 1]*u.pix)
    with pytest.raises(ValueError, match='At least one coordinate is not within the bounds of the map.'):
        sample_at_coords(aia171_test_map, point)


def test_contains_solar_center(aia171_test_map, all_off_disk_map, all_on_disk_map, straddles_limb_map, sub_smap):
    assert contains_solar_center(aia171_test_map)
    assert not contains_solar_center(all_off_disk_map)
    assert not contains_solar_center(all_on_disk_map)
    assert not contains_solar_center(straddles_limb_map)
    assert not contains_solar_center(sub_smap)


def test_verify_coordinate_helioprojective(aia171_test_map, all_off_disk_map, all_on_disk_map, straddles_limb_map,
                                           sub_smap, non_helioprojective_map, non_helioprojective_skycoord):
    # These should be helioprojective.
    _verify_coordinate_helioprojective(aia171_test_map.coordinate_frame)
    _verify_coordinate_helioprojective(all_off_disk_map.coordinate_frame)
    _verify_coordinate_helioprojective(all_on_disk_map.coordinate_frame)
    _verify_coordinate_helioprojective(straddles_limb_map.coordinate_frame)
    _verify_coordinate_helioprojective(sub_smap.coordinate_frame)
    # These are not.
    with pytest.raises(ValueError, match=r"HeliographicCarrington, .* Helioprojective"):
        _verify_coordinate_helioprojective(non_helioprojective_map.coordinate_frame)
    with pytest.raises(ValueError, match=r"ICRS, .* Helioprojective"):
        _verify_coordinate_helioprojective(non_helioprojective_skycoord)


def test_functions_raise_non_frame_coords(non_helioprojective_skycoord):
    with pytest.raises(ValueError, match=r"ICRS, .* Helioprojective"):
        solar_angular_radius(non_helioprojective_skycoord)
    with pytest.raises(ValueError, match=r"ICRS, .* Helioprojective"):
        coordinate_is_on_solar_disk(non_helioprojective_skycoord)


def test_functions_raise_non_frame_map(non_helioprojective_map):
    with pytest.raises(ValueError, match=r"HeliographicCarrington, .* Helioprojective"):
        contains_full_disk(non_helioprojective_map)
    with pytest.raises(ValueError, match=r"HeliographicCarrington, .* Helioprojective"):
        contains_solar_center(non_helioprojective_map)
    with pytest.raises(ValueError, match=r"HeliographicCarrington, .* Helioprojective"):
        is_all_off_disk(non_helioprojective_map)
    with pytest.raises(ValueError, match=r"HeliographicCarrington, .* Helioprojective"):
        contains_limb(non_helioprojective_map)
    with pytest.raises(ValueError, match=r"HeliographicCarrington, .* Helioprojective"):
        on_disk_bounding_coordinates(non_helioprojective_map)


def test_contains_coord(aia171_test_map):
    smap = aia171_test_map
    for coord in [smap.bottom_left_coord,
                  smap.top_right_coord,
                  SkyCoord(0*u.deg, 0*u.deg, frame=smap.coordinate_frame)]:
        assert contains_coordinate(smap, coord)

    assert not contains_coordinate(smap, SkyCoord(2000*u.arcsec, 2000*u.arcsec,
                                                  frame=smap.coordinate_frame))

    multi_coord = SkyCoord([0, 2000]*u.arcsec, [0, 2000]*u.arcsec,
                           frame=smap.coordinate_frame)
    assert (contains_coordinate(smap, multi_coord) == [True, False]).all()


@pytest.mark.parametrize(('x', 'y', 'sampled_x', 'sampled_y'),
                         [([1, 5], [1, 1], [1, 2, 3, 4, 5], [1, 1, 1, 1, 1]),
                          ([1, 5], [1, 2], [1, 2, 3, 3, 4, 5], [1, 1, 1, 2, 2, 2]),
                          ([1, 5], [1, 3], [1, 2, 2, 3, 4, 4, 5], [1, 1, 2, 2, 2, 3, 3]),
                          ([1, 5], [1.0, 4.0], [1, 2, 2, 3, 3, 4, 4, 5], [1, 1, 2, 2, 3, 3, 4, 4]),
                          ([1, 5], [1.0, 3.6], [1, 2, 2, 3, 3, 4, 5, 5], [1, 1, 2, 2, 3, 3, 3, 4]),
                          ([1, 5], [1.4, 3.6], [1, 1, 2, 3, 3, 4, 5, 5], [1, 2, 2, 2, 3, 3, 3, 4])])
def test_pixelate_coord_path(aia171_test_map, x, y, sampled_x, sampled_y):
    # Also test the x<->y transpose
    for xx, yy, sxx, syy in [(x, y, sampled_x, sampled_y), (y, x, sampled_y, sampled_x)]:
        # Using the AIA test map for a "real" WCS, but the actual WCS is irrelevant for this test
        line = aia171_test_map.wcs.pixel_to_world(xx, yy)
        sampled_coords = pixelate_coord_path(aia171_test_map, line)
        sampled_pixels = aia171_test_map.wcs.world_to_pixel(sampled_coords)
        assert np.allclose(sampled_pixels[0], sxx)
        assert np.allclose(sampled_pixels[1], syy)


@pytest.mark.parametrize(('x', 'y', 'sampled_x', 'sampled_y'),
                         [([1, 5], [1, 1], [1, 2, 3, 4, 5], [1, 1, 1, 1, 1]),
                          ([1, 5], [1, 2], [1, 2, 3, 4, 5], [1, 1, 1, 2, 2]),
                          ([1, 5], [1, 3], [1, 2, 3, 4, 5], [1, 1, 2, 2, 3]),
                          ([1, 5], [1.0, 4.0], [1, 2, 3, 4, 5], [1, 2, 2, 3, 4]),
                          ([1, 5], [1.0, 3.6], [1, 2, 3, 4, 5], [1, 2, 2, 3, 4]),
                          ([1, 5], [1.4, 3.6], [1, 2, 3, 4, 5], [1, 2, 2, 3, 4])])
def test_pixelate_coord_path_bresenham(aia171_test_map, x, y, sampled_x, sampled_y):
    # Also test the x<->y transpose
    for xx, yy, sxx, syy in [(x, y, sampled_x, sampled_y), (y, x, sampled_y, sampled_x)]:
        # Using the AIA test map for a "real" WCS, but the actual WCS is irrelevant for this test
        line = aia171_test_map.wcs.pixel_to_world(xx, yy)
        sampled_coords = pixelate_coord_path(aia171_test_map, line, bresenham=True)
        sampled_pixels = aia171_test_map.wcs.world_to_pixel(sampled_coords)
        assert np.allclose(sampled_pixels[0], sxx)
        assert np.allclose(sampled_pixels[1], syy)
