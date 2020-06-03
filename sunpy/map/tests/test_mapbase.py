"""
Test Generic Map
"""
import os
import tempfile
from unittest import mock

import matplotlib.pyplot as plt
import numpy as np
import pytest

import astropy.units as u
import astropy.wcs
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.io.fits.verify import VerifyWarning
from astropy.tests.helper import assert_quantity_allclose
from astropy.time import Time
from astropy.visualization import wcsaxes

import sunpy
import sunpy.coordinates
import sunpy.data.test
import sunpy.map
import sunpy.sun
from sunpy.coordinates import sun
from sunpy.time import parse_time
from sunpy.util import SunpyDeprecationWarning, SunpyUserWarning

testpath = sunpy.data.test.rootdir


@pytest.fixture
def hmi_test_map():
    return sunpy.map.Map(os.path.join(testpath, "resampled_hmi.fits"))


@pytest.fixture
def aia171_test_map():
    map = sunpy.map.Map(os.path.join(testpath, 'aia_171_level1.fits'))

    # Get rid of the blank keyword to prevent some astropy fits fixing warnings
    header = dict(map.meta)
    header.pop('blank')
    return sunpy.map.Map((map.data, header))


@pytest.fixture
def heliographic_test_map():
    map = sunpy.map.Map(os.path.join(testpath, 'heliographic_phase_map.fits.gz'))

    # Fix unit strings to prevent some astropy fits fixing warnings
    header = dict(map.meta)
    header['cunit1'] = 'deg'
    header['cunit2'] = 'deg'
    return sunpy.map.Map((map.data, header))


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
        'date-obs': '1970-01-01T00:00:00',
        'mjd-obs': 40587.0,
        'obsrvtry': 'Foo',
        'detector': 'bar',
        'wavelnth': 10,
        'waveunit': 'm'
    }
    return sunpy.map.Map((data, header))


@pytest.fixture
def simple_map():
    # A 3x3 map, with it's center at (0, 0), and scaled differently in
    # each direction
    data = np.arange(9).reshape((3, 3))
    ref_coord = SkyCoord(0.0, 0.0, frame='helioprojective', obstime='now', unit='deg')
    ref_pix = [1, 1] * u.pix
    scale = [2, 1] * u.arcsec / u.pix
    header = sunpy.map.make_fitswcs_header(data, ref_coord, reference_pixel=ref_pix, scale=scale)

    return sunpy.map.Map(data, header)


def test_fits_data_comparison(aia171_test_map):
    """Make sure the data is the same in pyfits and SunPy"""
    with pytest.warns(VerifyWarning, match="Invalid 'BLANK' keyword in header."):
        data = fits.open(os.path.join(testpath, 'aia_171_level1.fits'))[0].data
    np.testing.assert_allclose(aia171_test_map.data, data)


def test_get_item(generic_map):
    with pytest.raises(NotImplementedError):
        generic_map[10, 10]


def test_wcs(aia171_test_map):
    wcs = aia171_test_map.wcs
    assert isinstance(wcs, astropy.wcs.WCS)

    assert all(wcs.wcs.crpix - 1 ==
               [aia171_test_map.reference_pixel.x.value, aia171_test_map.reference_pixel.y.value])
    assert u.allclose(wcs.wcs.cdelt * (u.Unit(wcs.wcs.cunit[0])/u.pix),
                      u.Quantity(aia171_test_map.scale))
    assert u.allclose(wcs.wcs.crval * u.Unit(wcs.wcs.cunit[0]),
                      u.Quantity([aia171_test_map._reference_longitude, aia171_test_map._reference_latitude]))
    assert set(wcs.wcs.ctype) == {
        aia171_test_map.coordinate_system.axis1, aia171_test_map.coordinate_system.axis2}
    np.testing.assert_allclose(wcs.wcs.pc, aia171_test_map.rotation_matrix)


def test_header_immutability(aia171_test_map):
    # Check that accessing the wcs of a map doesn't modify the meta data
    assert 'KEYCOMMENTS' in aia171_test_map.meta
    aia171_test_map.wcs
    assert 'KEYCOMMENTS' in aia171_test_map.meta


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
    assert isinstance(generic_map.name, str)


def test_nickname(generic_map):
    assert generic_map.nickname == 'bar'


def test_nickname_set(generic_map):
    assert generic_map.nickname == 'bar'
    generic_map.nickname = 'hi'
    assert generic_map.nickname == 'hi'


def test_date(generic_map):
    assert isinstance(generic_map.date, Time)


def test_date_aia(aia171_test_map):
    assert aia171_test_map.date == parse_time('2011-02-15T00:00:00.34')


def test_detector(generic_map):
    assert generic_map.detector == 'bar'


def test_dsun(generic_map):
    with pytest.warns(SunpyUserWarning, match='Missing metadata for observer: assuming Earth-based observer.*'):
        assert_quantity_allclose(generic_map.dsun, sun.earth_distance(generic_map.date))


def test_rsun_meters(generic_map):
    assert generic_map.rsun_meters == sunpy.sun.constants.radius


def test_rsun_obs(generic_map):
    with pytest.warns(SunpyUserWarning, match='Missing metadata for solar radius'):
        assert generic_map.rsun_obs == sun.angular_radius(generic_map.date)


def test_coordinate_system(generic_map):
    assert generic_map.coordinate_system == ('HPLN-   ', 'HPLT-   ')


def test_carrington_longitude(generic_map):
    with pytest.warns(SunpyUserWarning, match='Missing metadata for observer: assuming Earth-based observer.*'):
        assert generic_map.carrington_longitude == sun.L0(generic_map.date)


def test_heliographic_latitude(generic_map):
    with pytest.warns(SunpyUserWarning, match='Missing metadata for observer: assuming Earth-based observer.*'):
        assert generic_map.heliographic_latitude == sun.B0(generic_map.date)


def test_heliographic_longitude(generic_map):
    with pytest.warns(SunpyUserWarning, match='Missing metadata for observer: assuming Earth-based observer.*'):
        assert generic_map.heliographic_longitude == 0.


def test_units(generic_map):
    generic_map.spatial_units == ('arcsec', 'arcsec')


def test_cmap(generic_map):
    assert generic_map.cmap == plt.get_cmap('gray')


def test_coordinate_frame(aia171_test_map):
    frame = aia171_test_map.coordinate_frame
    assert isinstance(frame, sunpy.coordinates.Helioprojective)
    assert frame.observer.lat == aia171_test_map.observer_coordinate.frame.lat
    assert frame.observer.lon == aia171_test_map.observer_coordinate.frame.lon
    assert frame.observer.radius == aia171_test_map.observer_coordinate.frame.radius
    assert frame.obstime == aia171_test_map.date


def test_heliographic_longitude_crln(hmi_test_map):
    assert_quantity_allclose(hmi_test_map.heliographic_longitude,
                             hmi_test_map.carrington_longitude - sun.L0(hmi_test_map.date),
                             rtol=1e-3)  # A tolerance is needed because L0 is for Earth, not SDO


def test_remove_observers(aia171_test_map):
    aia171_test_map._remove_existing_observer_location()
    with pytest.warns(SunpyUserWarning,
                      match='Missing metadata for observer: assuming Earth-based observer.*'):
        aia171_test_map.observer_coordinate


def test_partially_missing_observers(generic_map):
    generic_map.meta['hglt_obs'] = 0
    generic_map.meta['hgln_obs'] = 0
    generic_map.meta['crlt_obs'] = 0
    generic_map.meta['crln_obs'] = 0
    with pytest.warns(SunpyUserWarning,
                      match="Missing metadata for observer: assuming Earth-based observer.\n" +
                            "For frame 'heliographic_stonyhurst' the following metadata is missing: dsun_obs\n" +
                            "For frame 'heliographic_carrington' the following metadata is missing: dsun_obs\n"):
        generic_map.observer_coordinate

# ==============================================================================
# Test Rotation WCS conversion
# ==============================================================================


def test_rotation_matrix_pci_j(generic_map):
    np.testing.assert_allclose(generic_map.rotation_matrix, np.array([[0., -1.], [1., 0.]]))


def test_rotation_matrix_crota(aia171_test_map):
    np.testing.assert_allclose(aia171_test_map.rotation_matrix,
                               np.array([[9.99999943e-01, -3.38820761e-04],
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
        'NAXIS2': 6,
        'CUNIT1': 'arcsec',
        'CUNIT2': 'arcsec'
    }
    cd_map = sunpy.map.Map((data, header))
    np.testing.assert_allclose(cd_map.rotation_matrix, np.array([[0., -1.], [1., 0]]))


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
        'NAXIS2': 6,
        'CUNIT1': 'arcsec',
        'CUNIT2': 'arcsec'
    }
    cd_map = sunpy.map.Map((data, header))
    np.testing.assert_allclose(cd_map.rotation_matrix, np.array([[0., -1], [1., 0]]))


def test_swap_cd():
    amap = sunpy.map.Map(os.path.join(testpath, 'swap_lv1_20140606_000113.fits'))
    np.testing.assert_allclose(amap.rotation_matrix, np.array([[1., 0], [0, 1.]]))


def test_world_to_pixel(generic_map):
    """Make sure conversion from data units to pixels is internally
    consistent"""
    with pytest.warns(SunpyUserWarning, match='Missing metadata for observer'):
        # Note: FITS pixels start from 1,1
        test_pixel = generic_map.world_to_pixel(generic_map.reference_coordinate)
    assert_quantity_allclose(test_pixel, generic_map.reference_pixel)


def test_world_to_pixel_error(generic_map):
    strerr = 'world_to_pixel takes a Astropy coordinate frame or SkyCoord instance'
    with pytest.raises(ValueError, match=strerr):
        generic_map.world_to_pixel(1)


def test_save(aia171_test_map, generic_map):
    """Tests the map save function"""
    aiamap = aia171_test_map
    afilename = tempfile.NamedTemporaryFile(suffix='fits').name
    aiamap.save(afilename, filetype='fits', overwrite=True)
    loaded_save = sunpy.map.Map(afilename)
    assert isinstance(loaded_save, sunpy.map.sources.AIAMap)
    assert loaded_save.meta == aiamap.meta
    assert_quantity_allclose(loaded_save.data, aiamap.data)


def test_save_compressed(aia171_test_map, generic_map):
    """Tests the map save function"""
    aiamap = aia171_test_map
    afilename = tempfile.NamedTemporaryFile(suffix='fits').name
    aiamap.save(afilename, filetype='fits', hdu_type=fits.CompImageHDU, overwrite=True)
    loaded_save = sunpy.map.Map(afilename)
    # We expect that round tripping to CompImageHDU will change the header and
    # the data a little.
    assert isinstance(loaded_save, sunpy.map.sources.AIAMap)


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
        'NAXIS2': 6,
        'CUNIT1': 'arcsec',
        'CUNIT2': 'arcsec',
    }
    cd_map = sunpy.map.Map((data, header))
    assert cd_map.shifted_value[0].value == 0
    assert cd_map.shifted_value[1].value == 0


def test_shift_applied(generic_map):
    """Test that adding a shift actually updates the reference coordinate"""
    with pytest.warns(SunpyUserWarning, match='Missing metadata for observer'):
        original_reference_coord = (generic_map.reference_coordinate.Tx,
                                    generic_map.reference_coordinate.Ty)
    x_shift = 5 * u.arcsec
    y_shift = 13 * u.arcsec
    shifted_map = generic_map.shift(x_shift, y_shift)
    with pytest.warns(SunpyUserWarning, match='Missing metadata for observer'):
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


def test_corners(simple_map):
    # These are the centers of the corner pixels
    assert u.allclose(simple_map.top_right_coord.Tx, 2 * u.arcsec)
    assert u.allclose(simple_map.top_right_coord.Ty, 1 * u.arcsec)
    assert u.allclose(simple_map.bottom_left_coord.Tx, -2 * u.arcsec)
    assert u.allclose(simple_map.bottom_left_coord.Ty, -1 * u.arcsec)


def test_center(simple_map):
    assert u.allclose(simple_map.center.Tx, 0 * u.arcsec, atol=1e-26 * u.arcsec)
    assert u.allclose(simple_map.center.Ty, 0 * u.arcsec)


def test_dimensions(simple_map):
    assert simple_map.dimensions[0] == 3 * u.pix
    assert simple_map.dimensions[1] == 3 * u.pix


pixel_corners = [
    [([0, 0] * u.pix, [0, 0] * u.pix), np.array([[0]])],
    [([-1, -1] * u.pix, [0, 0] * u.pix), np.array([[0]])],
    # 0.5, 0.5 is the edge of the first pixel, so make sure
    # we don't include any other pixels
    [([0, 0] * u.pix, [0.5, 0.5] * u.pix), np.array([[0]])],
    [([0, 0] * u.pix, [0, 0.51] * u.pix), np.array([[0],
                                                    [3]])],
    [([0, 0] * u.pix, [0.51, 0] * u.pix), np.array([[0, 1]])],
    [([0, 0] * u.pix, [0.51, 0.51] * u.pix), np.array([[0, 1],
                                                       [3, 4]])],
    [([0.1, 0.1] * u.pix, [1.6, 1.4] * u.pix), np.array([[0, 1, 2],
                                                         [3, 4, 5]])],
    [([0, 0] * u.pix, [20, 20] * u.pix), np.array([[0, 1, 2],
                                                   [3, 4, 5],
                                                   [6, 7, 8]])],
]


@pytest.mark.parametrize('rect, submap_out', pixel_corners)
def test_submap_pixel(simple_map, rect, submap_out):
    # Check that result is the same specifying corners either way round
    for r in [dict(bottom_left=rect[0], top_right=rect[1]),
              dict(bottom_left=rect[1], top_right=rect[0])]:
        submap = simple_map.submap(**r)
        np.testing.assert_equal(submap.data, submap_out)


# The (0.5, 0.5) case is skipped as boundary points cannot reliably tested when
# converting to world coordinates due to round-off error when round-tripping
# through pixel_to_world -> world_to_pixel
@pytest.mark.parametrize('rect, submap_out', pixel_corners[:2] + pixel_corners[3:])
def test_submap_world(simple_map, rect, submap_out):
    # Check that coordinates behave the same way
    corner1 = simple_map.pixel_to_world(*rect[0])
    corner2 = simple_map.pixel_to_world(*rect[1])
    corners = simple_map.pixel_to_world(u.Quantity([rect[0][0], rect[1][0]]),
                                        u.Quantity([rect[0][1], rect[1][1]]))
    for r in [dict(bottom_left=corner1, top_right=corner2),
              dict(bottom_left=corner2, top_right=corner1),
              dict(bottom_left=corners, ),
              ]:
        submap = simple_map.submap(**r)
        np.testing.assert_equal(submap.data, submap_out)


# Check that submap works with units convertable to pix but that aren't pix
@pytest.mark.parametrize('unit', [u.pix, u.mpix * 1e3])
def test_submap_data_header(generic_map, unit):
    """Check data and header information for a submap"""
    width = generic_map.data.shape[1]
    height = generic_map.data.shape[0]

    # Create a submap of the top-right quadrant of the image
    submap = generic_map.submap([width / 2., height / 2.] * unit, top_right=[width, height] * unit)

    # Check to see if submap properties were updated properly
    assert submap.reference_pixel.x.value == generic_map.meta['crpix1'] - 1 - width / 2.
    assert submap.reference_pixel.y.value == generic_map.meta['crpix2'] - 1 - height / 2.
    assert submap.data.shape[1] == width / 2.
    assert submap.data.shape[0] == height / 2.

    # Check to see if header was updated
    assert submap.meta['naxis1'] == width / 2.
    assert submap.meta['naxis2'] == height / 2.

    # Check data
    assert (generic_map.data[height // 2:height, width // 2:width] == submap.data).all()


def test_reference_coordinate(simple_map):
    assert simple_map.reference_pixel.x == 1 * u.pix
    assert simple_map.reference_pixel.y == 1 * u.pix


@pytest.mark.parametrize('shape', [[1, 1], [6, 6]])
def test_resample(simple_map, shape):
    # Test resampling a 2x2 map
    resampled = simple_map.resample(shape * u.pix, method='linear')
    assert np.mean(resampled.data) == np.mean(simple_map.data)
    # Should be the mean of [0,1,2,3,4,5,6,7,8,9]
    if shape == [1, 1]:
        assert resampled.data == np.array([[4]])
    assert resampled.scale.axis1 == 3 / shape[0] * simple_map.scale.axis1
    assert resampled.scale.axis2 == 3 / shape[1] * simple_map.scale.axis2

    # Check that the corner coordinates of the input and output are the same
    resampled_lower_left = resampled.pixel_to_world(-0.5 * u.pix, -0.5 * u.pix)
    original_lower_left = simple_map.pixel_to_world(-0.5 * u.pix, -0.5 * u.pix)
    assert resampled_lower_left.Tx == original_lower_left.Tx
    assert resampled_lower_left.Ty == original_lower_left.Ty

    resampled_upper_left = resampled.pixel_to_world((shape[0] - 0.5) * u.pix,
                                                    (shape[1] - 0.5) * u.pix)
    original_upper_left = simple_map.pixel_to_world(2.5 * u.pix, 2.5 * u.pix)
    assert resampled_upper_left.Tx == original_upper_left.Tx
    assert resampled_upper_left.Ty == original_upper_left.Ty


resample_test_data = [('linear', (100, 200) * u.pixel), ('neighbor', (128, 256) * u.pixel),
                      ('nearest', (512, 128) * u.pixel), ('spline', (200, 200) * u.pixel)]


@pytest.mark.parametrize('sample_method, new_dimensions', resample_test_data)
def test_resample_dimensions(generic_map, sample_method, new_dimensions):
    """Check that resampled map has expected dimensions."""
    with pytest.warns(SunpyUserWarning, match='Missing metadata for observer'):
        resampled_map = generic_map.resample(new_dimensions, method=sample_method)
    assert resampled_map.dimensions[0] == new_dimensions[0]
    assert resampled_map.dimensions[1] == new_dimensions[1]


@pytest.mark.parametrize('sample_method, new_dimensions', resample_test_data)
def test_resample_metadata(generic_map, sample_method, new_dimensions):
    """
    Check that the resampled map has correctly adjusted metadata.
    """
    with pytest.warns(SunpyUserWarning, match='Missing metadata for observer'):
        resampled_map = generic_map.resample(new_dimensions, method=sample_method)
    assert float(resampled_map.meta['cdelt1']) / generic_map.meta['cdelt1'] \
        == float(generic_map.data.shape[1]) / resampled_map.data.shape[1]
    assert float(resampled_map.meta['cdelt2']) / generic_map.meta['cdelt2'] \
        == float(generic_map.data.shape[0]) / resampled_map.data.shape[0]
    assert resampled_map.meta['crpix1'] == (resampled_map.data.shape[1] + 1) / 2.
    assert resampled_map.meta['crpix2'] == (resampled_map.data.shape[0] + 1) / 2.
    with pytest.warns(SunpyUserWarning, match='Missing metadata for observer'):
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


def test_superpixel_err(generic_map):
    with pytest.raises(ValueError, match="Offset is strictly non-negative."):
        generic_map.superpixel((2, 2) * u.pix, offset=(-2, 2) * u.pix)


def calc_new_matrix(angle):
    c = np.cos(np.deg2rad(angle))
    s = np.sin(np.deg2rad(angle))
    return np.array([[c, -s], [s, c]])


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
    with pytest.warns(SunpyUserWarning, match='Missing metadata for observer'):
        rotated_map = generic_map.rotate(30*u.deg)
    # This tests that the reference pixel of the map is in the expected place.
    assert rotated_map.data.shape != generic_map.data.shape
    assert_quantity_allclose(u.Quantity(rotated_map.reference_pixel),
                             u.Quantity((5.04903811, 6.54903811), u.pix))


def test_rotate_recenter(generic_map):
    with pytest.warns(SunpyUserWarning, match='Missing metadata for observer'):
        rotated_map = generic_map.rotate(20 * u.deg, recenter=True)
    pixel_array_center = (np.flipud(rotated_map.data.shape) - 1) / 2.0

    assert_quantity_allclose(
        pixel_array_center * u.pix, u.Quantity(rotated_map.reference_pixel))


def test_rotate_crota_remove(aia171_test_map):
    rot_map = aia171_test_map.rotate()
    assert rot_map.meta.get('CROTA1', None) is None
    assert rot_map.meta.get('CROTA2', None) is None


def test_rotate_scale_cdelt(generic_map):
    with pytest.warns(SunpyUserWarning, match='Missing metadata for observer'):
        rot_map = generic_map.rotate(scale=10.)
    assert rot_map.meta['CDELT1'] == generic_map.meta['CDELT1'] / 10.
    assert rot_map.meta['CDELT2'] == generic_map.meta['CDELT2'] / 10.


def test_rotate_new_matrix(generic_map):
    # Rotate by CW90 to go from CCW 90 in generic map to CCW 180
    with pytest.warns(SunpyUserWarning, match='Missing metadata for observer'):
        rot_map = generic_map.rotate(rmatrix=np.array([[0, 1], [-1, 0]]))
    np.testing.assert_allclose(rot_map.rotation_matrix, np.array([[-1, 0], [0, -1]]))


def test_rotate_rmatrix_angle(generic_map):
    with pytest.raises(ValueError, match="You cannot specify both an angle and a rotation matrix."):
        generic_map.rotate(angle=5*u.deg, rmatrix=np.array([[1, 0], [0, 1]]))


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
    with pytest.warns(SunpyUserWarning, match='Missing metadata for observer'):
        out = generic_map.pixel_to_world(*u.Quantity(generic_map.reference_pixel))
    assert_quantity_allclose(out.Tx, 0*u.arcsec)
    assert_quantity_allclose(out.Ty, 0*u.arcsec)


def test_validate_meta(generic_map):
    """Check to see if_validate_meta displays an appropriate error"""
    with pytest.warns(SunpyUserWarning) as w:
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
        sunpy.map.Map((generic_map.data, bad_header))

    assert 'waveunit'.upper() in str(w[0].message)


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
    hdr['cunit1'] = 'arcsec'
    hdr['cunit2'] = 'arcsec'
    with pytest.warns(SunpyUserWarning, match='This file contains more than 2 dimensions.'):
        bad_map = sunpy.map.Map(bad_data, hdr)
    # Test fails if map.ndim > 2 and if the dimensions of the array are wrong.
    assert bad_map.ndim == 2
    assert_quantity_allclose(bad_map.dimensions, (5, 3) * u.pix)


def test_missing_metadata_warnings():
    # Checks that warnings for missing metadata are only raised once
    with pytest.warns(Warning) as record:
        header = {}
        header['cunit1'] = 'arcsec'
        header['cunit2'] = 'arcsec'
        array_map = sunpy.map.Map(np.random.rand(20, 15), header)
        array_map.peek()
    # There should be 2 warnings for missing metadata (obstime and observer location)
    assert len([w for w in record if w.category is SunpyUserWarning]) == 2


def test_fits_header(aia171_test_map):
    assert isinstance(aia171_test_map.fits_header, fits.Header)


def test_bad_coordframe_repr(generic_map):
    generic_map.meta['CTYPE1'] = "STUART1"
    generic_map.meta['CTYPE2'] = "STUART2"
    assert 'Unknown' in generic_map.__repr__()


def test_bad_header_final_fallback():
    # Checks that if a WCS cannot be constructed from the
    # header, a warning is raised and a simple WCS is created
    # instead
    header = {'cunit1': 'arcsec',
              'cunit2': 'arcsec',
              None: None,  # Cannot parse this into WCS
              }
    m = sunpy.map.GenericMap(np.zeros((10, 10)), header)
    with pytest.warns(UserWarning,
                      match="Unable to treat `.meta` as a FITS header, assuming a simple WCS."):
        m.wcs
    # Simple WCS validation
    assert list(m.wcs.wcs.ctype) == ['HPLN-', 'HPLT-']
    assert (m.wcs.wcs.crval == [0.0, 0.0]).all()
    assert (m.wcs.wcs.crpix == [5.5, 5.5]).all()
    assert (m.wcs.wcs.cdelt == [1.0, 1.0]).all()


def test_wcs_isot(aia171_test_map):
    # Check that a Map WCS returns the time as isot format
    assert aia171_test_map.wcs.to_header()['DATE-OBS'] == '2011-02-15T00:00:00.340'


def test_repr_html(aia171_test_map):
    html_string = aia171_test_map._repr_html_()
    assert isinstance(html_string, str)

    # Add a NaN value and check
    aia171_test_map.data[0, 0] = np.nan
    html_string = aia171_test_map._repr_html_()
    assert "Bad pixels are shown in red: 1 NaN" in html_string

    # Add a infinite value and check
    aia171_test_map.data[0, 0] = np.inf
    html_string = aia171_test_map._repr_html_()
    assert "Bad pixels are shown in red: 1 infinite" in html_string


def test_quicklook(aia171_test_map):
    with mock.patch('webbrowser.open_new_tab') as mockwbopen:
        aia171_test_map.quicklook()

    # Check that the mock web browser was opened with a file URL
    mockwbopen.assert_called_once()
    file_url = mockwbopen.call_args[0][0]
    assert file_url.startswith('file://')

    # Open the file specified in the URL and confirm that it contains the HTML
    with open(file_url[7:], 'r') as f:
        html_string = f.read()

        assert aia171_test_map._repr_html_() in html_string


@pytest.fixture
def generic_map2(generic_map):
    generic_map.meta["CTYPE1"] = "HPLN-TAN"
    generic_map.meta["CTYPE2"] = "HPLT-TAN"
    return generic_map


@pytest.fixture
def coords(generic_map2):
    bl_coord = SkyCoord(20, -10, unit=u.arcsec,
                        frame=generic_map2.coordinate_frame)

    tr_coord = SkyCoord(0, 10, unit=u.arcsec,
                        frame=generic_map2.coordinate_frame)

    bl_tr_coord = SkyCoord([20, 0], [-10, 10], unit=u.arcsec,
                           frame=generic_map2.coordinate_frame)

    return bl_coord, tr_coord, bl_tr_coord


bl_pix = [3, 2] * u.pix
tr_pix = [5, 4] * u.pix
width_pix = 2 * u.pix
height_pix = 2 * u.pix
width_deg = 20 * u.arcsec
height_deg = 20 * u.arcsec


def test_deprecated_submap_inputs(generic_map2, coords):
    bl_coord, tr_coord, bl_tr_coord = coords
    # deprecated
    with pytest.warns(SunpyDeprecationWarning):
        smap = generic_map2.submap(bl_coord, tr_coord)
    assert u.allclose(smap.dimensions, (3, 3) * u.pix)

    with pytest.warns(SunpyDeprecationWarning):
        smap = generic_map2.submap(bl_pix, tr_pix)
    assert u.allclose(smap.dimensions, (3, 3) * u.pix)

    # error
    with pytest.raises(TypeError, match="width must be specified as a keyword argument"):
        generic_map2.submap(bl_coord, width_deg, height_deg)

    with pytest.raises(TypeError, match="width must be specified as a keyword argument"):
        generic_map2.submap(bl_pix, width_pix, height_pix)

    with pytest.warns(SunpyDeprecationWarning):
        with pytest.raises(ValueError,
                           match="Either top_right alone or both width and height must be specified."):
            generic_map2.submap(bl_coord, width_deg, height=height_deg)

    with pytest.warns(SunpyDeprecationWarning):
        with pytest.raises(ValueError,
                           match="Either top_right alone or both width and height must be specified"):
            generic_map2.submap(bl_pix, width_pix, height=height_pix)

    with pytest.warns(SunpyDeprecationWarning):
        with pytest.raises(TypeError,
                           match="Invalid input, top_right must be of type SkyCoord or BaseCoordinateFrame."):
            generic_map2.submap(bl_coord, width_deg)

    with pytest.warns(SunpyDeprecationWarning):
        with pytest.raises(ValueError,
                           match="Either top_right alone or both width and height must be specified"):
            generic_map2.submap(bl_pix, width_pix, height=height_pix)


@pytest.mark.skip
def test_submap_kwarg_only_input_errors(generic_map2, coords):
    """
    This test replaces the one above when the deprecation period is over.
    """
    bl_coord, tr_coord, bl_tr_coord = coords
    inputs = (
        ((bl_coord, tr_coord), {}),
        ((bl_pix, tr_pix), {}),
        ((bl_coord, width_deg, height_deg), {}),
        ((bl_pix, width_pix, height_pix), {}),
        ((bl_coord, width_deg), {'height_deg': height_deg}),
        ((bl_pix, width_pix), {'height_pix': height_pix}),
    )

    for args, kwargs in inputs:
        with pytest.raises(TypeError, match="too many positional arguments"):
            generic_map2.submap(*args, **kwargs)


def test_submap_kwarg_only_input_errors(generic_map2, coords):
    bl_coord, tr_coord, bl_tr_coord = coords

    inputs = (
        dict(width=width_pix),
        dict(top_right=tr_pix, width=width_pix),
        dict(top_right=tr_pix, height=height_pix),
        dict(),  # Only post deprecation
    )
    for kwargs in inputs:
        with pytest.raises(ValueError, match="top_right alone or both width and height "
                                             "must be specified"):
            generic_map2.submap(bl_pix, **kwargs)

    with pytest.raises(TypeError, match="width and height must be a Quantity in units of pixels"):
        generic_map2.submap(bl_pix, width=width_deg, height=height_deg)

    with pytest.raises(TypeError,
                       match="top_right must be a Quantity in units of pixels."):
        generic_map2.submap([10, 10]*u.deg, top_right=[10, 10]*u.deg)

    with pytest.raises(ValueError, match=r"must have shape \(2\, \)"):
        generic_map2.submap(10*u.pix, top_right=10*u.pix)

    with pytest.raises(ValueError, match=r"must have shape \(2\, \)"):
        generic_map2.submap([10, 10, 10]*u.pix, top_right=[10, 10, 10]*u.pix)

    with pytest.raises(u.UnitsError):
        generic_map2.submap([10, 10]*u.deg, width=10*u.km, height=10*u.J)

    with pytest.raises(ValueError,
                       match="either bottom_left and top_right or bottom_left and height and width should be provided"):
        generic_map2.submap(SkyCoord([10, 10, 10]*u.deg, [10, 10, 10]*u.deg,
                                     frame=generic_map2.coordinate_frame))


def test_submap_inputs(generic_map2, coords):
    bl_coord, tr_coord, bl_tr_coord = coords

    inputs = (
        ((bl_coord,), dict(top_right=tr_coord)),
        ((bl_coord,), dict(width=width_deg, height=height_deg)),
        ((bl_tr_coord,), {}),
        ((bl_pix,), dict(top_right=tr_pix)),
        ((bl_pix,), dict(width=width_pix, height=height_pix)),
    )

    for args, kwargs in inputs:
        smap = generic_map2.submap(*args, **kwargs)
        assert u.allclose(smap.dimensions, (3, 3) * u.pix)
