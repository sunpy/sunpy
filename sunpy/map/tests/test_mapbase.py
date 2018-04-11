# -*- coding: utf-8 -*-
"""
Test Generic Map
"""
from __future__ import absolute_import

import os
import pytest
import datetime
import warnings
import tempfile

import numpy as np

import astropy.wcs
from astropy.io import fits
import astropy.units as u
from astropy.tests.helper import assert_quantity_allclose
from astropy.visualization import wcsaxes
from astropy.coordinates import SkyCoord
import matplotlib.pyplot as plt

import sunpy
import sunpy.sun
import sunpy.map
import sunpy.coordinates
import sunpy.data.test
from sunpy.time import parse_time

from sunpy.extern import six

testpath = sunpy.data.test.rootdir


@pytest.fixture
def aia171_test_map():
    return sunpy.map.Map(os.path.join(testpath, 'aia_171_level1.fits'))


@pytest.fixture
def heliographic_test_map():
    return sunpy.map.Map(os.path.join(testpath, 'heliographic_phase_map.fits.gz'))


@pytest.fixture
def aia171_test_map_with_mask(aia171_test_map):
    shape = aia171_test_map.data.shape
    mask = np.zeros_like(aia171_test_map.data, dtype=bool)
    mask[0:shape[0] // 2, 0:shape[1] // 2] = True
    return sunpy.map.Map(np.ma.array(aia171_test_map.data, mask=mask), aia171_test_map.meta)


@pytest.fixture
def generic_map():
    data = np.ones([6, 6], dtype=np.float64)
    header = {
        'CRVAL1': 0,
        'CRVAL2': 0,
        'CRPIX1': 5,
        'CRPIX2': 5,
        'CDELT1': 10,
        'CDELT2': 10,
        'CUNIT1': 'arcsec',
        'CUNIT2': 'arcsec',
        'PC1_1': 0,
        'PC1_2': -1,
        'PC2_1': 1,
        'PC2_2': 0,
        'NAXIS1': 6,
        'NAXIS2': 6,
        'date-obs': '1970/01/01T00:00:00',
        'obsrvtry': 'Foo',
        'detector': 'bar',
        'wavelnth': 10,
        'waveunit': 'm'
    }
    return sunpy.map.Map((data, header))


def test_fits_data_comparison(aia171_test_map):
    """Make sure the data is the same in pyfits and SunPy"""
    data = fits.open(os.path.join(testpath, 'aia_171_level1.fits'))[0].data
    np.testing.assert_allclose(aia171_test_map.data, data)


def test_get_item(generic_map):
    with pytest.raises(NotImplementedError):
        generic_map[10, 10]


def test_wcs(aia171_test_map):
    wcs = aia171_test_map.wcs
    assert isinstance(wcs, astropy.wcs.WCS)

    assert all(wcs.wcs.crpix ==
               [aia171_test_map.reference_pixel.x.value, aia171_test_map.reference_pixel.y.value])
    assert all(wcs.wcs.cdelt == [aia171_test_map.scale.axis1.value,
                                 aia171_test_map.scale.axis1.value])
    assert all(
        wcs.wcs.crval ==
        [aia171_test_map._reference_longitude.value, aia171_test_map._reference_latitude.value])
    assert set(wcs.wcs.ctype) == set(
        [aia171_test_map.coordinate_system.axis1, aia171_test_map.coordinate_system.axis2])
    np.testing.assert_allclose(wcs.wcs.pc, aia171_test_map.rotation_matrix)
    assert set(wcs.wcs.cunit) == set([u.Unit(a) for a in aia171_test_map.spatial_units])


def test_dtype(generic_map):
    assert generic_map.dtype == np.float64


def test_size(generic_map):
    assert generic_map.size == 36 * u.pix


def test_min(generic_map):
    assert generic_map.min() == 1


def test_max(generic_map):
    assert generic_map.max() == 1


def test_mean(generic_map):
    assert generic_map.mean() == 1


def test_std(generic_map):
    assert generic_map.std() == 0


# ==============================================================================
# Test the default value of a load of properties
# TODO: Test the header keyword extraction
# ==============================================================================
def test_name(generic_map):
    assert isinstance(generic_map.name, six.string_types)


def test_nickname(generic_map):
    assert generic_map.nickname == 'bar'


def test_nickname_set(generic_map):
    assert generic_map.nickname == 'bar'
    generic_map.nickname = 'hi'
    assert generic_map.nickname == 'hi'


def test_date(generic_map):
    assert isinstance(generic_map.date, datetime.datetime)


def test_date_aia(aia171_test_map):
    assert aia171_test_map.date == parse_time('2011-02-15T00:00:00.34')


def test_detector(generic_map):
    assert generic_map.detector == 'bar'


def test_dsun(generic_map):
    assert generic_map.dsun == sunpy.coordinates.get_sunearth_distance(generic_map.date).to(u.m)


def test_rsun_meters(generic_map):
    assert generic_map.rsun_meters == sunpy.sun.constants.radius


def test_rsun_obs(generic_map):
    assert generic_map.rsun_obs == sunpy.sun.solar_semidiameter_angular_size(generic_map.date)


def test_coordinate_system(generic_map):
    assert generic_map.coordinate_system == ('HPLN-   ', 'HPLT-   ')


def test_carrington_longitude(generic_map):
    assert generic_map.carrington_longitude == sunpy.coordinates.get_sun_L0(generic_map.date)


def test_heliographic_latitude(generic_map):
    assert generic_map.heliographic_latitude == sunpy.coordinates.get_sun_B0(generic_map.date)


def test_heliographic_longitude(generic_map):
    assert generic_map.heliographic_longitude == 0.


def test_units(generic_map):
    generic_map.spatial_units == ('arcsec', 'arcsec')


def test_coordinate_frame(aia171_test_map):
    frame = aia171_test_map.coordinate_frame
    assert isinstance(frame, sunpy.coordinates.Helioprojective)
    assert frame.observer.lat == aia171_test_map.observer_coordinate.frame.lat
    assert frame.observer.lon == aia171_test_map.observer_coordinate.frame.lon
    assert frame.observer.radius == aia171_test_map.observer_coordinate.frame.radius
    assert frame.obstime == aia171_test_map.date


# ==============================================================================
# Test Rotation WCS conversion
# ==============================================================================
def test_rotation_matrix_pci_j(generic_map):
    np.testing.assert_allclose(generic_map.rotation_matrix, np.matrix([[0., -1.], [1., 0.]]))


def test_rotation_matrix_crota(aia171_test_map):
    np.testing.assert_allclose(aia171_test_map.rotation_matrix,
                               np.matrix([[9.99999943e-01, -3.38820761e-04],
                                          [3.38820761e-04, 9.99999943e-01]]))


def test_rotation_matrix_cd_cdelt():
    data = np.ones([6, 6], dtype=np.float64)
    header = {
        'CRVAL1': 0,
        'CRVAL2': 0,
        'CRPIX1': 5,
        'CRPIX2': 5,
        'CDELT1': 10,
        'CDELT2': 9,
        'CD1_1': 0,
        'CD1_2': -9,
        'CD2_1': 10,
        'CD2_2': 0,
        'NAXIS1': 6,
        'NAXIS2': 6
    }
    cd_map = sunpy.map.Map((data, header))
    np.testing.assert_allclose(cd_map.rotation_matrix, np.matrix([[0., -1.], [1., 0]]))


def test_rotation_matrix_cd_cdelt_square():
    data = np.ones([6, 6], dtype=np.float64)
    header = {
        'CRVAL1': 0,
        'CRVAL2': 0,
        'CRPIX1': 5,
        'CRPIX2': 5,
        'CDELT1': 10,
        'CDELT2': 10,
        'CD1_1': 0,
        'CD1_2': -10,
        'CD2_1': 10,
        'CD2_2': 0,
        'NAXIS1': 6,
        'NAXIS2': 6
    }
    cd_map = sunpy.map.Map((data, header))
    np.testing.assert_allclose(cd_map.rotation_matrix, np.matrix([[0., -1], [1., 0]]))


def test_swap_cd():
    amap = sunpy.map.Map(os.path.join(testpath, 'swap_lv1_20140606_000113.fits'))
    np.testing.assert_allclose(amap.rotation_matrix, np.matrix([[1., 0], [0, 1.]]))


def test_data_range(generic_map):
    """Make sure xrange and yrange work"""
    assert_quantity_allclose(
        (generic_map.xrange[1] - generic_map.xrange[0]).to(u.arcsec).value,
        generic_map.meta['cdelt1'] * generic_map.meta['naxis1']
    )
    assert_quantity_allclose(
        (generic_map.yrange[1] - generic_map.yrange[0]).to(u.arcsec).value,
        generic_map.meta['cdelt2'] * generic_map.meta['naxis2']
    )

    # the weird unit-de-unit thing here is to work around and inconsistency in
    # the way np.average works with astropy 1.3 and 2.0dev
    assert_quantity_allclose(np.average(u.Quantity(generic_map.xrange).value) * u.arcsec,
                             generic_map.center.Tx)
    assert_quantity_allclose(np.average(u.Quantity(generic_map.yrange).value) * u.arcsec,
                             generic_map.center.Ty)


def test_world_to_pixel(generic_map):
    """Make sure conversion from data units to pixels is internally
    consistent"""
    # Note: FITS pixels start from 1,1
    test_pixel = generic_map.world_to_pixel(generic_map.reference_coordinate, origin=1)
    assert_quantity_allclose(test_pixel, generic_map.reference_pixel)


def test_save(generic_map):
    """Tests the map save function"""
    aiamap = aia171_test_map()
    afilename = tempfile.NamedTemporaryFile(suffix='fits').name
    aiamap.save(afilename, filetype='fits', clobber=True)
    loaded_save = sunpy.map.Map(afilename)
    assert isinstance(loaded_save, sunpy.map.sources.AIAMap)
    assert loaded_save.meta == aiamap.meta
    assert_quantity_allclose(loaded_save.data, aiamap.data)


def test_default_shift():
    """Test that the default shift is zero"""
    data = np.ones([6, 6], dtype=np.float64)
    header = {
        'CRVAL1': 0,
        'CRVAL2': 0,
        'CRPIX1': 5,
        'CRPIX2': 5,
        'CDELT1': 10,
        'CDELT2': 9,
        'CD1_1': 0,
        'CD1_2': -9,
        'CD2_1': 10,
        'CD2_2': 0,
        'NAXIS1': 6,
        'NAXIS2': 6
    }
    cd_map = sunpy.map.Map((data, header))
    assert cd_map.shifted_value[0].value == 0
    assert cd_map.shifted_value[1].value == 0


def test_shift_applied(generic_map):
    """Test that adding a shift actually updates the reference coordinate"""
    original_reference_coord = (generic_map.reference_coordinate.Tx,
                                generic_map.reference_coordinate.Ty)
    x_shift = 5 * u.arcsec
    y_shift = 13 * u.arcsec
    shifted_map = generic_map.shift(x_shift, y_shift)
    assert shifted_map.reference_coordinate.Tx - x_shift == original_reference_coord[0]
    assert shifted_map.reference_coordinate.Ty - y_shift == original_reference_coord[1]
    crval1 = ((generic_map.meta.get('crval1') * generic_map.spatial_units[0] +
               shifted_map.shifted_value[0]).to(shifted_map.spatial_units[0])).value
    assert shifted_map.meta.get('crval1') == crval1
    crval2 = ((generic_map.meta.get('crval2') * generic_map.spatial_units[1] +
               shifted_map.shifted_value[1]).to(shifted_map.spatial_units[1])).value
    assert shifted_map.meta.get('crval2') == crval2


def test_set_shift(generic_map):
    """Test that previously applied shift is stored in the shifted_value property"""
    x_shift = 5 * u.arcsec
    y_shift = 13 * u.arcsec
    shifted_map = generic_map.shift(x_shift, y_shift)
    resultant_shift = shifted_map.shifted_value
    assert resultant_shift[0] == x_shift
    assert resultant_shift[1] == y_shift


def test_shift_history(generic_map):
    """Test the shifted_value is added to a non-zero previous shift"""
    x_shift1 = 5 * u.arcsec
    y_shift1 = 13 * u.arcsec
    shifted_map1 = generic_map.shift(x_shift1, y_shift1)

    x_shift2 = -28.5 * u.arcsec
    y_shift2 = 120 * u.arcsec
    final_shifted_map = shifted_map1.shift(x_shift2, y_shift2)

    resultant_shift = final_shifted_map.shifted_value
    assert resultant_shift[0] == x_shift1 + x_shift2
    assert resultant_shift[1] == y_shift1 + y_shift2


def test_submap(generic_map):
    """Check data and header information for a submap"""
    width = generic_map.data.shape[1]
    height = generic_map.data.shape[0]

    # Create a submap of the top-right quadrant of the image
    submap = generic_map.submap([width / 2., height / 2.] * u.pix, [width, height] * u.pix)

    # Check to see if submap properties were updated properly
    assert submap.reference_pixel.x.value == generic_map.meta['crpix1'] - width / 2.
    assert submap.reference_pixel.y.value == generic_map.meta['crpix2'] - height / 2.
    assert submap.data.shape[1] == width / 2.
    assert submap.data.shape[0] == height / 2.

    # Check to see if header was updated
    assert submap.meta['naxis1'] == width / 2.
    assert submap.meta['naxis2'] == height / 2.

    # Check data
    assert (generic_map.data[height // 2:height, width // 2:width] == submap.data).all()


resample_test_data = [('linear', (100, 200) * u.pixel), ('neighbor', (128, 256) * u.pixel),
                      ('nearest', (512, 128) * u.pixel), ('spline', (200, 200) * u.pixel)]


@pytest.mark.parametrize('sample_method, new_dimensions', resample_test_data)
def test_resample_dimensions(generic_map, sample_method, new_dimensions):
    """Check that resampled map has expected dimensions."""
    resampled_map = generic_map.resample(new_dimensions, method=sample_method)
    assert resampled_map.dimensions[0] == new_dimensions[0]
    assert resampled_map.dimensions[1] == new_dimensions[1]


@pytest.mark.parametrize('sample_method, new_dimensions', resample_test_data)
def test_resample_metadata(generic_map, sample_method, new_dimensions):
    """
    Check that the resampled map has correctly adjusted metadata.
    """
    resampled_map = generic_map.resample(new_dimensions, method=sample_method)
    assert float(resampled_map.meta['cdelt1']) / generic_map.meta['cdelt1'] \
        == float(generic_map.data.shape[1]) / resampled_map.data.shape[1]
    assert float(resampled_map.meta['cdelt2']) / generic_map.meta['cdelt2'] \
        == float(generic_map.data.shape[0]) / resampled_map.data.shape[0]
    assert resampled_map.meta['crpix1'] == (resampled_map.data.shape[1] + 1) / 2.
    assert resampled_map.meta['crpix2'] == (resampled_map.data.shape[0] + 1) / 2.
    assert resampled_map.meta['crval1'] == generic_map.center.Tx.value
    assert resampled_map.meta['crval2'] == generic_map.center.Ty.value
    for key in generic_map.meta:
        if key not in ('cdelt1', 'cdelt2', 'crpix1', 'crpix2', 'crval1', 'crval2'):
            assert resampled_map.meta[key] == generic_map.meta[key]


def test_superpixel(aia171_test_map, aia171_test_map_with_mask):
    dimensions = (2, 2) * u.pix
    superpixel_map_sum = aia171_test_map.superpixel(dimensions)
    assert_quantity_allclose(superpixel_map_sum.dimensions[1],
                             aia171_test_map.dimensions[1] / dimensions[1] * u.pix)
    assert_quantity_allclose(superpixel_map_sum.dimensions[0],
                             aia171_test_map.dimensions[0] / dimensions[0] * u.pix)
    assert_quantity_allclose(superpixel_map_sum.data[0][0],
                             (aia171_test_map.data[0][0] + aia171_test_map.data[0][1] +
                              aia171_test_map.data[1][0] + aia171_test_map.data[1][1]))

    superpixel_map_avg = aia171_test_map.superpixel(dimensions, func=np.mean)
    assert_quantity_allclose(superpixel_map_avg.dimensions[1],
                             aia171_test_map.dimensions[1] / dimensions[1] * u.pix)
    assert_quantity_allclose(superpixel_map_avg.dimensions[0],
                             aia171_test_map.dimensions[0] / dimensions[0] * u.pix)
    assert_quantity_allclose(superpixel_map_avg.data[0][0],
                             (aia171_test_map.data[0][0] + aia171_test_map.data[0][1] +
                              aia171_test_map.data[1][0] + aia171_test_map.data[1][1]) / 4.0)

    # Test that the mask is respected
    superpixel_map_sum = aia171_test_map_with_mask.superpixel(dimensions)
    assert superpixel_map_sum.mask is not None
    assert_quantity_allclose(superpixel_map_sum.mask.shape[0],
                             aia171_test_map.dimensions[1] / dimensions[1])
    assert_quantity_allclose(superpixel_map_sum.mask.shape[1],
                             aia171_test_map.dimensions[0] / dimensions[0])

    # Test that the offset is respected
    superpixel_map_sum = aia171_test_map_with_mask.superpixel(dimensions, offset=(1, 1) * u.pix)
    assert_quantity_allclose(superpixel_map_sum.dimensions[1],
                             aia171_test_map.dimensions[1] / dimensions[1] * u.pix - 1 * u.pix)
    assert_quantity_allclose(superpixel_map_sum.dimensions[0],
                             aia171_test_map.dimensions[0] / dimensions[0] * u.pix - 1 * u.pix)

    dimensions = (7, 9) * u.pix
    superpixel_map_sum = aia171_test_map_with_mask.superpixel(dimensions, offset=(4, 4) * u.pix)
    assert_quantity_allclose(
        superpixel_map_sum.dimensions[0],
        np.int((aia171_test_map.dimensions[0] / dimensions[0]).value) * u.pix - 1 * u.pix)
    assert_quantity_allclose(
        superpixel_map_sum.dimensions[1],
        np.int((aia171_test_map.dimensions[1] / dimensions[1]).value) * u.pix - 1 * u.pix)


def calc_new_matrix(angle):
    c = np.cos(np.deg2rad(angle))
    s = np.sin(np.deg2rad(angle))
    return np.matrix([[c, -s], [s, c]])


def test_rotate(aia171_test_map):
    rotated_map_1 = aia171_test_map.rotate(20 * u.deg)
    rotated_map_2 = rotated_map_1.rotate(20 * u.deg)
    np.testing.assert_allclose(rotated_map_1.rotation_matrix,
                               np.dot(aia171_test_map.rotation_matrix, calc_new_matrix(20).T))
    np.testing.assert_allclose(rotated_map_2.rotation_matrix,
                               np.dot(aia171_test_map.rotation_matrix, calc_new_matrix(40).T))

    # Rotation of a map by a non-integral multiple of 90 degrees expands the map
    # and assigns the value of 0 to corner pixels. This results in a reduction
    # of the mean for a map of all non-negative values.
    assert rotated_map_2.data.shape > rotated_map_1.data.shape > aia171_test_map.data.shape
    np.testing.assert_allclose(rotated_map_1.data[0, 0], 0., atol=1e-7)
    np.testing.assert_allclose(rotated_map_2.data[0, 0], 0., atol=1e-7)
    assert rotated_map_2.mean() < rotated_map_1.mean() < aia171_test_map.mean()

    rotated_map_3 = aia171_test_map.rotate(0 * u.deg, scale=1.5)
    assert rotated_map_3.mean() > aia171_test_map.mean()

    # Mean and std should be equal when angle of rotation is integral multiple
    # of 90 degrees for a square map
    rotated_map_4 = aia171_test_map.rotate(90 * u.deg, scale=1.5)
    np.testing.assert_allclose(rotated_map_3.mean(), rotated_map_4.mean(), rtol=1e-3)
    np.testing.assert_allclose(rotated_map_3.std(), rotated_map_4.std(), rtol=1e-3)
    rotated_map_5 = aia171_test_map.rotate(180 * u.deg, scale=1.5)
    np.testing.assert_allclose(rotated_map_3.mean(), rotated_map_5.mean(), rtol=1e-3)
    np.testing.assert_allclose(rotated_map_3.std(), rotated_map_5.std(), rtol=2e-3)

    # Rotation of a rectangular map by a large enough angle will change which dimension is larger
    aia171_test_map_crop = aia171_test_map.submap(
        SkyCoord(
            [[0, 0], [1000, 400]] * u.arcsec, frame=aia171_test_map.coordinate_frame))

    aia171_test_map_crop_rot = aia171_test_map_crop.rotate(60 * u.deg)
    assert aia171_test_map_crop.data.shape[0] < aia171_test_map_crop.data.shape[1]
    assert aia171_test_map_crop_rot.data.shape[0] > aia171_test_map_crop_rot.data.shape[1]

    # Same test as above, to test the other direction
    aia171_test_map_crop = aia171_test_map.submap(
        SkyCoord(
            [[0, 0], [400, 1000]] * u.arcsec, frame=aia171_test_map.coordinate_frame))
    aia171_test_map_crop_rot = aia171_test_map_crop.rotate(60 * u.deg)
    assert aia171_test_map_crop.data.shape[0] > aia171_test_map_crop.data.shape[1]
    assert aia171_test_map_crop_rot.data.shape[0] < aia171_test_map_crop_rot.data.shape[1]


def test_rotate_pad_crpix(generic_map):
    rotated_map = generic_map.rotate(30*u.deg)
    # This tests that the reference pixel of the map is in the expected place.
    assert rotated_map.data.shape != generic_map.data.shape
    assert_quantity_allclose(u.Quantity(rotated_map.reference_pixel),
                             u.Quantity((6.049038105675565, 7.5490381056760265), u.pix))


def test_rotate_recenter(generic_map):
    rotated_map = generic_map.rotate(20 * u.deg, recenter=True)
    pixel_array_center = (np.flipud(rotated_map.data.shape) - 1) / 2.0

    assert_quantity_allclose(
        (pixel_array_center + 1) * u.pix,  # FITS indexes from 1
        u.Quantity(rotated_map.reference_pixel))


def test_rotate_crota_remove(aia171_test_map):
    rot_map = aia171_test_map.rotate()
    assert rot_map.meta.get('CROTA1', None) is None
    assert rot_map.meta.get('CROTA2', None) is None


def test_rotate_scale_cdelt(generic_map):
    rot_map = generic_map.rotate(scale=10.)
    assert rot_map.meta['CDELT1'] == generic_map.meta['CDELT1'] / 10.
    assert rot_map.meta['CDELT2'] == generic_map.meta['CDELT2'] / 10.


def test_rotate_new_matrix(generic_map):
    # Rotate by CW90 to go from CCW 90 in generic map to CCW 180
    rot_map = generic_map.rotate(rmatrix=np.matrix([[0, 1], [-1, 0]]))
    np.testing.assert_allclose(rot_map.rotation_matrix, np.matrix([[-1, 0], [0, -1]]))


def test_rotate_rmatrix_angle(generic_map):
    with pytest.raises(ValueError):
        generic_map.rotate(angle=5, rmatrix=np.matrix([[1, 0], [0, 1]]))


def test_rotate_invalid_order(generic_map):
    with pytest.raises(ValueError):
        generic_map.rotate(order=6)
    with pytest.raises(ValueError):
        generic_map.rotate(order=-1)


def test_as_mpl_axes_aia171(aia171_test_map):
    ax = plt.subplot(projection=aia171_test_map)
    assert isinstance(ax, wcsaxes.WCSAxes)
    # This test doesn't work, it seems that WCSAxes copies or changes the WCS
    # object.
    #  assert ax.wcs is aia171_test_map.wcs
    assert all([ct1 == ct2 for ct1, ct2 in zip(ax.wcs.wcs.ctype, aia171_test_map.wcs.wcs.ctype)])
    # Map adds these attributes, so we use them to check.
    assert hasattr(ax.wcs, 'heliographic_observer')


def test_pixel_to_world_no_projection(generic_map):
    out = generic_map.pixel_to_world(*u.Quantity(generic_map.reference_pixel)+1*u.pix, origin=1)
    assert_quantity_allclose(out.Tx, -10*u.arcsec)
    assert_quantity_allclose(out.Ty, 10*u.arcsec)


def test_validate_meta(generic_map):
    """Check to see if_validate_meta displays an appropriate error"""
    with warnings.catch_warnings(record=True) as w:
        bad_header = {
            'CRVAL1': 0,
            'CRVAL2': 0,
            'CRPIX1': 5,
            'CRPIX2': 5,
            'CDELT1': 10,
            'CDELT2': 10,
            'CUNIT1': 'ARCSEC',
            'CUNIT2': 'ARCSEC',
            'PC1_1': 0,
            'PC1_2': -1,
            'PC2_1': 1,
            'PC2_2': 0,
            'NAXIS1': 6,
            'NAXIS2': 6,
            'date-obs': '1970/01/01T00:00:00',
            'obsrvtry': 'Foo',
            'detector': 'bar',
            'wavelnth': 10,
            'waveunit': 'ANGSTROM'
        }
        bad_map = sunpy.map.Map((generic_map.data, bad_header))
        for count, meta_property in enumerate(('cunit1', 'cunit2', 'waveunit')):
            assert meta_property.upper() in str(w[count].message)


# Heliographic Map Tests


def test_hg_coord(heliographic_test_map):
    assert heliographic_test_map.coordinate_system[0] == "CRLN-CAR"
    assert heliographic_test_map.coordinate_system[1] == "CRLT-CAR"
    assert isinstance(heliographic_test_map.coordinate_frame,
                      sunpy.coordinates.HeliographicCarrington)


def test_hg_pix_to_data(heliographic_test_map):
    out = heliographic_test_map.pixel_to_world(180 * u.pix, 90 * u.pix)
    assert isinstance(out, SkyCoord)
    assert isinstance(out.frame, sunpy.coordinates.HeliographicCarrington)
    assert_quantity_allclose(out.lon, 0 * u.deg)
    assert_quantity_allclose(out.lat, 0 * u.deg)


def test_hg_data_to_pix(heliographic_test_map):
    out = heliographic_test_map.world_to_pixel(
        SkyCoord(
            0 * u.deg, 0 * u.deg, frame=heliographic_test_map.coordinate_frame))
    assert_quantity_allclose(out[0], 180 * u.pix)
    assert_quantity_allclose(out[1], 90 * u.pix)

# Heliocentric Map Tests


def test_hc_warn():
    data = np.ones([6, 6], dtype=np.float64)
    header = {
        'CRVAL1': 0,
        'CRVAL2': 0,
        'CRPIX1': 5,
        'CRPIX2': 5,
        'CDELT1': 10,
        'CDELT2': 10,
        'CUNIT1': 'km',
        'CUNIT2': 'km',
        'CTYPE1': 'SOLX    ',
        'CTYPE2': 'SOLY    ',
        'PC1_1': 0,
        'PC1_2': -1,
        'PC2_1': 1,
        'PC2_2': 0,
        'NAXIS1': 6,
        'NAXIS2': 6,
        'date-obs': '1970/01/01T00:00:00',
        'obsrvtry': 'Foo',
        'detector': 'bar',
        'wavelnth': 10,
        'waveunit': 'm'
    }

    with pytest.warns(UserWarning):
        sunpy.map.Map((data, header))

# Dimension testing


def test_more_than_two_dimensions():
    """Checks to see if an appropriate error is raised when a FITS with more than two dimensions is
    loaded.  We need to load a >2-dim dataset with a TELESCOP header"""

    # Data crudely represnts 4 stokes, 4 wavelengths with Y,X of 3 and 5.
    bad_data = np.random.rand(4, 4, 3, 5)
    hdr = fits.Header()
    hdr['TELESCOP'] = 'XXX'
    bad_map = sunpy.map.Map(bad_data, hdr)
    # Test fails if map.ndim > 2 and if the dimensions of the array are wrong.
    assert bad_map.ndim is 2
    assert_quantity_allclose(bad_map.dimensions, (5, 3) * u.pix)


def test_missing_metadata_warnings():
    # Checks that warnings for missing metadata are only raised once
    with pytest.warns(Warning) as record:
        array_map = sunpy.map.Map(np.random.rand(20, 15), {})
        array_map.peek()
    # There should be 4 warnings for missing metadata
    assert len([w for w in record if 'Missing metadata' in str(w)]) == 4
