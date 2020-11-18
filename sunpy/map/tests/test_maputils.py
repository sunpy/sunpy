import os

import numpy as np
import pytest

import astropy.units as u
from astropy.coordinates import BaseCoordinateFrame, SkyCoord

import sunpy.data.test
import sunpy.map
from sunpy.coordinates import HeliographicStonyhurst
from sunpy.coordinates.frames import HeliographicCarrington
from sunpy.coordinates.utils import GreatArc
from sunpy.map.maputils import (
    _verify_coordinate_helioprojective,
    all_coordinates_from_map,
    all_pixel_indices_from_map,
    contains_full_disk,
    contains_limb,
    contains_solar_center,
    coordinate_is_on_solar_disk,
    is_all_off_disk,
    is_all_on_disk,
    map_edges,
    on_disk_bounding_coordinates,
    sample_at_coords,
    solar_angular_radius,
)

testpath = sunpy.data.test.rootdir


@pytest.fixture
def aia171_test_map():
    return sunpy.map.Map(os.path.join(testpath, 'aia_171_level1.fits'))


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


def test_map_edges(all_off_disk_map):
    edges = map_edges(all_off_disk_map)
    assert type(edges) is tuple
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


def test_data_at_coordinates(aia171_test_map, aia_test_arc):
    data = sample_at_coords(aia171_test_map, aia_test_arc.coordinates())
    pixels = np.asarray(np.rint(
        aia171_test_map.world_to_pixel(aia_test_arc.coordinates())), dtype=int)
    x = pixels[0, :]
    y = pixels[1, :]
    intensity_along_arc = aia171_test_map.data[y, x]
    np.testing.assert_almost_equal(data[0], intensity_along_arc[0], decimal=1)
    np.testing.assert_almost_equal(data[-1], intensity_along_arc[-1], decimal=1)


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
