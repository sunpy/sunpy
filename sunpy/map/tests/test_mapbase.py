# -*- coding: utf-8 -*-
"""
Test Generic Map
"""
from __future__ import absolute_import

import os

import numpy as np

import pytest

import astropy.wcs
from astropy.io import fits
import astropy.units as u

import sunpy
import sunpy.sun
import sunpy.map
import sunpy.data.test

testpath = sunpy.data.test.rootdir

@pytest.fixture
def aia171_test_map():
    return sunpy.map.Map(os.path.join(testpath, 'aia_171_level1.fits'))


@pytest.fixture
def generic_map():
    data = np.ones([6,6], dtype=np.float64)
    header = {'CRVAL1': 0,
              'CRVAL2': 0,
              'CRPIX1': 5,
              'CRPIX2': 5,
              'CDELT1': 10,
              'CDELT2': 10,
              'PC1_1': 0,
              'PC1_2': -1,
              'PC2_1': 1,
              'PC2_2': 0,
              'NAXIS1': 6,
              'NAXIS2': 6}
    return sunpy.map.Map((data, header))


#@pytest.fixture
#def aia171_test_map_large():
#    return sunpy.map.Map(sunpy.AIA_171_IMAGE)


def test_fits_data_comparison(aia171_test_map):
    """Make sure the data is the same in pyfits and SunPy"""
    data = fits.open(os.path.join(testpath, 'aia_171_level1.fits'))[0].data
    np.testing.assert_allclose(aia171_test_map.data, data)


def test_get_item(generic_map):
    with pytest.raises(NotImplementedError):
        generic_map[10,10]


def test_repr_no_obs(generic_map):
    assert generic_map.__repr__() == 'array([[ 1.,  1.,  1.,  1.,  1.,  1.],\n       [ 1.,  1.,  1.,  1.,  1.,  1.],\n       [ 1.,  1.,  1.,  1.,  1.,  1.],\n       [ 1.,  1.,  1.,  1.,  1.,  1.],\n       [ 1.,  1.,  1.,  1.,  1.,  1.],\n       [ 1.,  1.,  1.,  1.,  1.,  1.]])'


def test_repr_obs(aia171_test_map):
    assert aia171_test_map.__repr__() == 'SunPy AIAMap\n---------\nObservatory:\t SDO\nInstrument:\t AIA_3\nDetector:\t AIA\nMeasurement:\t 171\nObs Date:\t 2011-02-15T00:00:00.34\ndt:\t\t 2.000191\nDimension:\t [128, 128]\n[dx, dy] =\t [19.183648, 19.183648]\n\narray([[-1.25,  0.  ,  1.  , ...,  0.  ,  0.5 , -0.75],\n       [ 0.75, -0.25, -0.5 , ...,  0.25,  0.  , -0.25],\n       [ 0.  ,  0.5 ,  1.75, ...,  0.  ,  0.5 ,  0.  ],\n       ..., \n       [ 1.  ,  0.25, -0.25, ...,  0.  ,  0.  ,  0.  ],\n       [-0.25,  0.  , -0.5 , ...,  0.75, -0.75,  0.  ],\n       [ 0.75,  1.5 , -0.75, ...,  0.  , -0.5 ,  0.5 ]])'

def test_wcs(aia171_test_map):
    wcs = aia171_test_map.wcs
    assert isinstance(wcs, astropy.wcs.WCS)

    assert all(wcs.wcs.crpix == [aia171_test_map.reference_pixel['x'],
                                 aia171_test_map.reference_pixel['y']])
    assert all(wcs.wcs.cdelt == [aia171_test_map.scale['x'],
                                 aia171_test_map.scale['y']])
    assert all(wcs.wcs.crval == [aia171_test_map.reference_coordinate['x'],
                                 aia171_test_map.reference_coordinate['y']])
    assert set(wcs.wcs.ctype) == set([aia171_test_map.coordinate_system['x'],
                             aia171_test_map.coordinate_system['y']])
    np.testing.assert_allclose(wcs.wcs.pc, aia171_test_map.rotation_matrix)
    assert set(wcs.wcs.cunit) == set([u.Unit(a) for a in [aia171_test_map.units['x'],
                                                          aia171_test_map.units['y']]])

def test_dtype(generic_map):
    assert generic_map.dtype == np.float64


def test_size(generic_map):
    assert generic_map.size == 36


def test_min(generic_map):
    assert generic_map.min() == 1


def test_max(generic_map):
    assert generic_map.max() == 1


def test_mean(generic_map):
    assert generic_map.mean() == 1


def test_std(generic_map):
    assert generic_map.std() == 0


#==============================================================================
# Test the default value of a load of properties
# TODO: Test the header keyword extraction
#==============================================================================
def test_name(generic_map):
    assert generic_map.name == ' '


def test_name_set(generic_map):
    assert generic_map.name == ' '
    generic_map.name = 'hi'
    assert generic_map.name == 'hi'


def test_nickname(generic_map):
    assert generic_map.nickname == ''


def test_nickname_set(generic_map):
    assert generic_map.nickname == ''
    generic_map.nickname = 'hi'
    assert generic_map.nickname == 'hi'


def test_date(generic_map):
    assert generic_map.date is 'now'


def test_date_aia(aia171_test_map):
    assert aia171_test_map.date == '2011-02-15T00:00:00.34'


def test_detector(generic_map):
    assert generic_map.detector == ''


def test_dsun(generic_map):
    assert generic_map.dsun == sunpy.sun.sunearth_distance(generic_map.date).to(u.m)


def test_rsun_meters(generic_map):
    assert generic_map.rsun_meters == sunpy.sun.constants.radius


def test_rsun_arcseconds(generic_map):
    assert generic_map.rsun_arcseconds == sunpy.sun.solar_semidiameter_angular_size(generic_map.date).value


def test_coordinate_system(generic_map):
    assert generic_map.coordinate_system == {'x':'HPLN-TAN', 'y': 'HPLT-TAN'}


def test_carrington_longitude(generic_map):
    assert generic_map.carrington_longitude == (sunpy.sun.heliographic_solar_center(generic_map.date))[0]


def test_heliographic_latitude(generic_map):
    assert generic_map.heliographic_latitude == (sunpy.sun.heliographic_solar_center(generic_map.date))[1]


def test_heliographic_longitude(generic_map):
    assert generic_map.heliographic_longitude == 0.


def test_units(generic_map):
    generic_map.units == {'x': 'arcsec', 'y': 'arcsec'}


#==============================================================================
# Test Rotation WCS conversion
#==============================================================================
def test_rotation_matrix_pci_j(generic_map):
    np.testing.assert_allclose (generic_map.rotation_matrix,
                                np.matrix([[0., -1.], [1., 0.]]))


def test_rotation_matrix_crota(aia171_test_map):
    np.testing.assert_allclose (aia171_test_map.rotation_matrix,
                                np.matrix([[9.99999943e-01, -3.38820761e-04],
                                           [3.38820761e-04, 9.99999943e-01]]))


def test_rotation_matrix_cd_cdelt():
    data = np.ones([6,6], dtype=np.float64)
    header = {'CRVAL1': 0,
              'CRVAL2': 0,
              'CRPIX1': 5,
              'CRPIX2': 5,
              'CDELT1': 10,
              'CDELT2': 9,
              'CD1_1': 0,
              'CD1_2': -10,
              'CD2_1': 9,
              'CD2_2': 0,
              'NAXIS1': 6,
              'NAXIS2': 6}
    cd_map = sunpy.map.Map((data, header))
    np.testing.assert_allclose(cd_map.rotation_matrix, np.matrix([[0., -1.], [1., 0]]))

def test_rotation_matrix_cd_cdelt_square():
    data = np.ones([6,6], dtype=np.float64)
    header = {'CRVAL1': 0,
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
              'NAXIS2': 6}
    cd_map = sunpy.map.Map((data, header))
    np.testing.assert_allclose(cd_map.rotation_matrix, np.matrix([[0., -1], [1., 0]]))

def test_swap_cd():
    amap = sunpy.map.Map(os.path.join(testpath, 'swap_lv1_20140606_000113.fits'))
    np.testing.assert_allclose(amap.rotation_matrix, np.matrix([[1., 0], [0, 1.]]))


def test_data_range(generic_map):
    """Make sure xrange and yrange work"""
    assert generic_map.xrange[1] - generic_map.xrange[0] == generic_map.meta['cdelt1'] * generic_map.meta['naxis1']
    assert generic_map.yrange[1] - generic_map.yrange[0] == generic_map.meta['cdelt2'] * generic_map.meta['naxis2']

    assert np.average(generic_map.xrange) == generic_map.center['x']
    assert np.average(generic_map.yrange) == generic_map.center['y']


def test_data_to_pixel(generic_map):
    """Make sure conversion from data units to pixels is accurate"""
    # Check conversion of reference pixel
    # Note: FITS pixels starts from 1,1
    assert generic_map.data_to_pixel(generic_map.meta['crval1'], 'x') == generic_map.meta['crpix1'] - 1
    assert generic_map.data_to_pixel(generic_map.meta['crval2'], 'y') == generic_map.meta['crpix2'] - 1

    # Check conversion of map center
    assert generic_map.data_to_pixel(generic_map.center['x'], 'x') == (generic_map.meta['naxis1'] - 1) / 2.
    assert generic_map.data_to_pixel(generic_map.center['y'], 'y') == (generic_map.meta['naxis2'] - 1) / 2.

    # Check conversion of map edges
    # Note: data coords are at pixel centers, so edges are 0.5 pixels wider
    assert generic_map.data_to_pixel(generic_map.xrange[0], 'x') == 0. - 0.5
    assert generic_map.data_to_pixel(generic_map.yrange[0], 'y') == 0. - 0.5
    assert generic_map.data_to_pixel(generic_map.xrange[1], 'x') == (generic_map.meta['naxis1'] - 1) + 0.5
    assert generic_map.data_to_pixel(generic_map.yrange[1], 'y') == (generic_map.meta['naxis2'] - 1) + 0.5


def test_submap(generic_map):
    """Check data and header information for a submap"""
    width = generic_map.shape[1]
    height = generic_map.shape[0]

    # Create a submap of the top-right quadrant of the image
    submap = generic_map.submap([height/2.,height], [width/2.,width],
                                units='pixels')

    # Expected offset for center
    offset = {
        "x": generic_map.meta['crpix1'] - width / 2.,
        "y": generic_map.meta['crpix2'] - height / 2.,
    }

    # Check to see if submap properties were updated properly
    assert submap.reference_pixel['x'] == offset['x']
    assert submap.reference_pixel['y'] == offset['y']
    assert submap.shape[0] == width / 2.
    assert submap.shape[1] == height / 2.

    # Check to see if header was updated
    assert submap.meta['naxis1'] == width / 2.
    assert submap.meta['naxis2'] == height / 2.

    # Check data
    assert (generic_map.data[height/2:height,
                             width/2:width] == submap.data).all()


resample_test_data = [('linear', (100, 200)),
                      ('neighbor', (128, 256)),
                      ('nearest', (512, 128)),
                      ('spline', (200, 200))]

@pytest.mark.parametrize('sample_method, new_dimensions', resample_test_data)
def test_resample_dimensions(generic_map, sample_method, new_dimensions):
    """Check that resampled map has expected dimensions."""
    resampled_map = generic_map.resample(new_dimensions, method=sample_method)
    assert resampled_map.shape[1] == new_dimensions[0]
    assert resampled_map.shape[0] == new_dimensions[1]


@pytest.mark.parametrize('sample_method, new_dimensions', resample_test_data)
def test_resample_metadata(generic_map, sample_method, new_dimensions):
    """
    Check that the resampled map has correctly adjusted metadata.
    """
    resampled_map = generic_map.resample(new_dimensions, method=sample_method)
    assert float(resampled_map.meta['cdelt1']) / generic_map.meta['cdelt1'] \
        == float(generic_map.shape[1]) / resampled_map.shape[1]
    assert float(resampled_map.meta['cdelt2']) / generic_map.meta['cdelt2'] \
        == float(generic_map.shape[0]) / resampled_map.shape[0]
    assert resampled_map.meta['crpix1'] == (resampled_map.shape[1] + 1) / 2.
    assert resampled_map.meta['crpix2'] == (resampled_map.shape[0] + 1) / 2.
    assert resampled_map.meta['crval1'] == generic_map.center['x']
    assert resampled_map.meta['crval2'] == generic_map.center['y']
    for key in generic_map.meta:
        if key not in ('cdelt1', 'cdelt2', 'crpix1', 'crpix2',
                       'crval1', 'crval2'):
            assert resampled_map.meta[key] == generic_map.meta[key]


def test_superpixel(aia171_test_map):
    dimensions = (2, 2)
    superpixel_map_sum = aia171_test_map.superpixel(dimensions)
    assert superpixel_map_sum.shape[0] == aia171_test_map.shape[0]/dimensions[1]
    assert superpixel_map_sum.shape[1] == aia171_test_map.shape[1]/dimensions[0]
    assert superpixel_map_sum.data[0][0] == (aia171_test_map.data[0][0] +
                                             aia171_test_map.data[0][1] +
                                             aia171_test_map.data[1][0] +
                                             aia171_test_map.data[1][1])

    dimensions = (2, 2)
    superpixel_map_avg = aia171_test_map.superpixel(dimensions, 'average')
    assert superpixel_map_avg.shape[0] == aia171_test_map.shape[0]/dimensions[1]
    assert superpixel_map_avg.shape[1] == aia171_test_map.shape[1]/dimensions[0]
    assert superpixel_map_avg.data[0][0] == (aia171_test_map.data[0][0] +
                                             aia171_test_map.data[0][1] +
                                             aia171_test_map.data[1][0] +
                                             aia171_test_map.data[1][1])/4.0


def calc_new_matrix(angle):
    c = np.cos(np.deg2rad(angle))
    s = np.sin(np.deg2rad(angle))
    return np.matrix([[c, -s], [s, c]])


def test_rotate(aia171_test_map):
    rotated_map_1 = aia171_test_map.rotate(20)
    rotated_map_2 = rotated_map_1.rotate(20)
    np.testing.assert_allclose(rotated_map_1.rotation_matrix,
                               np.dot(aia171_test_map.rotation_matrix,
                                      calc_new_matrix(20).T))
    np.testing.assert_allclose(rotated_map_2.rotation_matrix,
                               np.dot(aia171_test_map.rotation_matrix,
                                      calc_new_matrix(40).T))

    # Rotation of a map by a non-integral multiple of 90 degrees expands the map
    # and assigns the value of 0 to corner pixels. This results in a reduction
    # of the mean for a map of all non-negative values.
    assert rotated_map_2.shape > rotated_map_1.shape > aia171_test_map.shape
    np.testing.assert_allclose(rotated_map_1.data[0, 0], 0., atol=1e-7)
    np.testing.assert_allclose(rotated_map_2.data[0, 0], 0., atol=1e-7)
    assert rotated_map_2.mean() < rotated_map_1.mean() < aia171_test_map.mean()

    rotated_map_3 = aia171_test_map.rotate(0, scale=1.5)
    assert rotated_map_3.mean() > aia171_test_map.mean()

    # Mean and std should be equal when angle of rotation is integral multiple
    # of 90 degrees for a square map
    rotated_map_4 = aia171_test_map.rotate(90, scale=1.5)
    np.testing.assert_allclose(rotated_map_3.mean(), rotated_map_4.mean(), rtol=1e-3)
    np.testing.assert_allclose(rotated_map_3.std(), rotated_map_4.std(), rtol=1e-3)
    rotated_map_5 = aia171_test_map.rotate(180, scale=1.5)
    np.testing.assert_allclose(rotated_map_3.mean(), rotated_map_5.mean(), rtol=1e-3)
    np.testing.assert_allclose(rotated_map_3.std(), rotated_map_5.std(), rtol=2e-3)

    # Rotation of a rectangular map by a large enough angle will change which dimension is larger
    aia171_test_map_crop = aia171_test_map.submap([0, 1000], [0, 400])
    aia171_test_map_crop_rot = aia171_test_map_crop.rotate(60)
    assert aia171_test_map_crop.shape[0] < aia171_test_map_crop.shape[1]
    assert aia171_test_map_crop_rot.shape[0] > aia171_test_map_crop_rot.shape[1]

    # Same test as above, to test the other direction
    aia171_test_map_crop = aia171_test_map.submap([0, 400], [0, 1000])
    aia171_test_map_crop_rot = aia171_test_map_crop.rotate(60)
    assert aia171_test_map_crop.shape[0] > aia171_test_map_crop.shape[1]
    assert aia171_test_map_crop_rot.shape[0] < aia171_test_map_crop_rot.shape[1]


def test_rotate_recenter(aia171_test_map):
    array_center = (np.array(aia171_test_map.data.shape)-1)/2.0

    # New image center in data coordinates
    new_center = np.asarray((200, 100))

    rotated_map = aia171_test_map.rotate(20, rotation_center=new_center, recenter=True)

    # Retrieve pixel coordinates for the centers from the new map
    new_x = rotated_map.data_to_pixel(new_center[0], 'x')
    new_y = rotated_map.data_to_pixel(new_center[1], 'y')

    # The new desired image center should be in the map center
    new_array_center = (np.array(rotated_map.data.shape)-1)/2.0
    np.testing.assert_allclose((new_y, new_x), new_array_center)


def test_rotate_crota_remove(aia171_test_map):
    rot_map = aia171_test_map.rotate()
    assert rot_map.meta.get('CROTA1', None) is None
    assert rot_map.meta.get('CROTA2', None) is None


def test_rotate_scale_cdelt(generic_map):
    rot_map = generic_map.rotate(scale=10.)
    assert rot_map.meta['CDELT1'] == generic_map.meta['CDELT1']/10.
    assert rot_map.meta['CDELT2'] == generic_map.meta['CDELT2']/10.


def test_rotate_new_matrix(generic_map):
    # Rotate by CW90 to go from CCW 90 in generic map to CCW 180
    rot_map = generic_map.rotate(rmatrix=np.matrix([[0, 1], [-1, 0]]))
    np.testing.assert_allclose(rot_map.rotation_matrix, np.matrix([[-1, 0], [0, -1]]))


def test_rotate_rmatrix_angle(generic_map):
    with pytest.raises(ValueError):
        generic_map.rotate(angle=5, rmatrix=np.matrix([[1,0], [0, 1]]))


def test_rotate_invalid_order(generic_map):
    with pytest.raises(ValueError):
        generic_map.rotate(order=6)
    with pytest.raises(ValueError):
        generic_map.rotate(order=-1)
