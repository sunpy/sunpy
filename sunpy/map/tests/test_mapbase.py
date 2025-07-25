"""
Test Generic Map
"""
import re
import tempfile
from copy import deepcopy

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pytest
from hypothesis import HealthCheck, given, settings
from matplotlib.figure import Figure
from matplotlib.transforms import Affine2D

import astropy.units as u
import astropy.wcs
from astropy.coordinates import Latitude, SkyCoord
from astropy.io import fits
from astropy.io.fits.verify import VerifyWarning
from astropy.tests.helper import assert_quantity_allclose
from astropy.visualization import wcsaxes
from astropy.wcs import InconsistentAxisTypesError

import sunpy
import sunpy.coordinates
import sunpy.map
import sunpy.sun
from sunpy.coordinates import HeliographicCarrington, HeliographicStonyhurst, sun
from sunpy.data.test import get_dummy_map_from_header, get_test_filepath
from sunpy.image.resample import reshape_image_to_4d_superpixel
from sunpy.image.transform import _rotation_registry
from sunpy.map.mapbase import GenericMap
from sunpy.map.sources import AIAMap
from sunpy.tests.helpers import asdf_entry_points, figure_test
from sunpy.time import parse_time
from sunpy.util import SunpyUserWarning
from sunpy.util.exceptions import SunpyDeprecationWarning, SunpyMetadataWarning
from sunpy.util.metadata import ModifiedItem
from .strategies import matrix_meta


def test_fits_data_comparison(aia171_test_map):
    """Make sure the data is the same when read with astropy.io.fits and sunpy"""
    with pytest.warns(VerifyWarning, match="Invalid 'BLANK' keyword in header."):
        hdulist = fits.open(get_test_filepath('aia_171_level1.fits'))
    np.testing.assert_allclose(aia171_test_map.data, hdulist[0].data)
    hdulist.close()


def test_header_fits_io():
    with pytest.warns(VerifyWarning, match="Invalid 'BLANK' keyword in header."):
        with fits.open(get_test_filepath('aia_171_level1.fits')) as hdu:
            AIAMap(hdu[0].data, hdu[0].header)


def test_get_item(generic_map):
    with pytest.raises(NotImplementedError):
        generic_map[10, 10]


def test_wcs(aia171_test_map):
    wcs = aia171_test_map.wcs
    assert isinstance(wcs, astropy.wcs.WCS)

    assert wcs.array_shape == aia171_test_map.data.shape

    assert all(wcs.wcs.crpix - 1 ==
               [aia171_test_map.reference_pixel.x.value, aia171_test_map.reference_pixel.y.value])
    assert u.allclose(wcs.wcs.cdelt * (u.Unit(wcs.wcs.cunit[0])/u.pix),
                      u.Quantity(aia171_test_map.scale))
    assert u.allclose(wcs.wcs.crval * u.Unit(wcs.wcs.cunit[0]),
                      u.Quantity([aia171_test_map._reference_longitude, aia171_test_map._reference_latitude]))
    assert set(wcs.wcs.ctype) == {
        aia171_test_map.coordinate_system.axis1, aia171_test_map.coordinate_system.axis2}
    np.testing.assert_allclose(wcs.wcs.pc, aia171_test_map.rotation_matrix)


def test_wcs_pv():
    # Test that PVi_m values are preserved in the reconstructed WCS
    zpn_header = {
        'ctype1': 'HPLN-ZPN',
        'ctype2': 'HPLT-ZPN',
        'cunit1': 'arcsec',
        'cunit2': 'arcsec',
        'pv1_0': 0,
        'pv1_1': 0,
        'pv1_2': 90,
        'pv1_3': 180,
        'pv2_1': 1,
        'pv2_5': 0.2,
        'pv2_10': 0.1,
        'date-obs': '2025-01-01',
        'hglt_obs': 0,
        'hgln_obs': 0,
        'dsun_obs': 1e13,
    }
    zpn_map = sunpy.map.Map((np.zeros((10, 10)), zpn_header))
    pv_values = zpn_map.wcs.wcs.get_pv()
    assert len(pv_values) == 7
    assert pv_values[0] == (1, 0, 0)
    assert pv_values[1] == (1, 1, 0)
    assert pv_values[2] == (1, 2, 90)
    assert pv_values[3] == (1, 3, 180)
    assert pv_values[4] == (2, 1, 1.0)
    assert pv_values[5] == (2, 5, 0.2)
    assert pv_values[6] == (2, 10, 0.1)


def test_wcs_cache(aia171_test_map):
    wcs1 = aia171_test_map.wcs
    wcs2 = aia171_test_map.wcs
    # Check that without any changes to the header, retrieving the wcs twice
    # returns the same object instead of recomputing the wcs
    assert wcs1 is wcs2

    # Change the header and make sure the wcs is re-computed
    new_crpix = 20
    assert new_crpix != wcs2.wcs.crpix[0]
    aia171_test_map.meta['crpix1'] = new_crpix

    new_wcs = aia171_test_map.wcs
    assert new_wcs.wcs.crpix[0] == new_crpix


def test_wcs_error_not_cached(aia171_test_map):
    # Create a cached value for the property
    _ = aia171_test_map.wcs

    # Modify the WCS in a bad way
    aia171_test_map.meta['ctype1'] = 'HPLN-ARC'

    # Try and fail to recalculate the property
    with pytest.raises(InconsistentAxisTypesError):
        _ = aia171_test_map.wcs

    # Try again and fail again to recalculate the property
    with pytest.raises(InconsistentAxisTypesError):
        _ = aia171_test_map.wcs


def test_obs_coord_cache(aia171_test_map):
    coord1 = aia171_test_map.observer_coordinate
    coord2 = aia171_test_map.observer_coordinate
    assert coord1 is coord2

    # Change metadata, and check that the coordinate changes
    aia171_test_map.meta['haex_obs'] += 10
    new_coord = aia171_test_map.observer_coordinate
    assert new_coord.lon != coord2.lon
    assert new_coord.lat != coord2.lat
    assert new_coord.radius != coord2.radius


def test_header_immutability(aia171_test_map):
    # Check that accessing the wcs of a map doesn't modify the meta data
    assert 'KEYCOMMENTS' in aia171_test_map.meta
    aia171_test_map.wcs
    assert 'KEYCOMMENTS' in aia171_test_map.meta


def test_dtype(generic_map):
    assert generic_map.dtype == np.float64


def test_min(generic_map):
    assert generic_map.min() == 0


def test_max(generic_map):
    assert generic_map.max() == 35


def test_mean(generic_map):
    assert generic_map.mean() == 17.5


def test_std(generic_map):
    np.testing.assert_allclose(generic_map.std(), 10.388294694831615)


def test_unit(generic_map):
    assert generic_map.unit == u.DN / u.s
    generic_map.meta['bunit'] = 'not a unit'
    with pytest.warns(SunpyMetadataWarning, match='Could not parse unit string "not a unit"'):
        assert generic_map.unit is None


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


date_dict = {'DATE-AVG': parse_time('2020-01-01'),
             'DATE-OBS': parse_time('2020-02-01'),
             'DATE-BEG': parse_time('2020-03-01'),
             'DATE-END': parse_time('2020-03-03')}
date_begend = date_dict['DATE-BEG'] + (date_dict['DATE-END'] - date_dict['DATE-BEG']) / 2


@pytest.mark.parametrize(("keys", "expected_date"),
                         [(['DATE-AVG', 'DATE-OBS', 'DATE-BEG', 'DATE-END'], date_dict['DATE-OBS']),
                          (['DATE-AVG', 'DATE-BEG', 'DATE-END'], date_dict['DATE-BEG']),
                          (['DATE-BEG', 'DATE-END'], date_dict['DATE-BEG']),
                          (['DATE-BEG'], date_dict['DATE-BEG']),
                          (['DATE-END'], date_dict['DATE-END']),
                          ([], 'now')
                          ])
def test_date(generic_map, keys, expected_date):
    # Remove pre-existing date keys
    for key in date_dict:
        generic_map.meta.pop(key, None)
    # Add new date keys
    for key in keys:
        generic_map.meta[key] = date_dict[key].isot
    # Check date is the correct value
    if expected_date == 'now':
        expected_date = parse_time('now')
        # Check equal to within a tolerance as parse_time('now') is run
        # at slightly different times in .date and the line above
        with pytest.warns(SunpyMetadataWarning, match='Missing metadata for observation time'):
            assert generic_map.date - expected_date < 1*u.s
    else:
        assert generic_map.date == expected_date


def test_date_scale(generic_map):
    # Check that default time scale is UTC
    assert 'timesys' not in generic_map.meta
    assert generic_map.date.scale == 'utc'
    generic_map.meta['timesys'] = 'tai'
    assert generic_map.date.scale == 'tai'


def test_date_aia(aia171_test_map):
    assert aia171_test_map.date == parse_time('2011-02-15T00:00:00.34')


def test_detector(generic_map):
    assert generic_map.detector == 'bar'


def test_timeunit(generic_map):
    assert generic_map.timeunit == u.Unit('s')
    generic_map.meta['timeunit'] = 'h'
    assert generic_map.timeunit == u.Unit('h')


def test_exposure_time(generic_map):
    exptime = 2 * u.s
    generic_map.meta['exptime'] = exptime.to_value('s')
    assert generic_map.exposure_time == exptime
    exptime = 3 * u.s
    # XPOSURE should take priority over EXPTIME
    generic_map.meta['xposure'] = exptime.to_value('s')
    assert generic_map.exposure_time == exptime
    del generic_map.meta['exptime']
    del generic_map.meta['xposure']
    assert generic_map.exposure_time is None
    # Test that an exposure time of 0.0 s does not yield None
    generic_map.meta['exptime'] = 0.0
    assert generic_map.exposure_time == 0.0 * u.s


def test_dsun(generic_map):
    assert_quantity_allclose(generic_map.dsun, sun.earth_distance(generic_map.date))


def test_rsun_meters(generic_map):
    assert generic_map.rsun_meters == sunpy.sun.constants.radius


def test_rsun_obs_without_rsun_ref(generic_map):
    assert_quantity_allclose(generic_map.rsun_obs,
                             sun.angular_radius(generic_map.date))


def test_rsun_obs_with_rsun_ref(generic_map):
    generic_map.meta['rsun_ref'] = sunpy.sun.constants.radius.to_value(u.m)
    # The following should not raise a warning because we can calculate it exactly
    assert_quantity_allclose(generic_map.rsun_obs, sun.angular_radius(generic_map.date))


def test_coordinate_system(generic_map):
    assert generic_map.coordinate_system == ('HPLN-TAN', 'HPLT-TAN')


def test_default_coordinate_system(generic_map):
    generic_map.meta.pop('ctype1')
    with pytest.warns(SunpyMetadataWarning, match='Missing CTYPE1 from metadata'):
        assert generic_map.coordinate_system == ('HPLN-TAN', 'HPLT-TAN')

    generic_map.meta.pop('ctype2')
    generic_map.meta['ctype1'] = 'HPLN-TAN'
    with pytest.warns(SunpyMetadataWarning, match='Missing CTYPE2 from metadata'):
        assert generic_map.coordinate_system == ('HPLN-TAN', 'HPLT-TAN')


@pytest.mark.skipif(pytest.__version__ < "8.0.0", reason="pytest >= 8.0.0 raises two warnings for this test")
def test_coordinate_system_solar_x_solar_y(generic_map):
    generic_map.meta['ctype1'] = 'SOLAR-X'
    generic_map.meta['ctype2'] = 'SOLAR-Y'
    with pytest.warns(SunpyDeprecationWarning, match="CTYPE1 value 'solar-x'/'solar_x' is deprecated") :
        with pytest.warns(SunpyDeprecationWarning, match="CTYPE2 value 'solar-y'/'solar_y' is deprecated") :
            assert generic_map.coordinate_system == ('HPLN-TAN', 'HPLT-TAN')


def test_carrington_longitude(generic_map):
    assert u.allclose(generic_map.carrington_longitude, sun.L0(generic_map.date))


def test_heliographic_latitude(generic_map):
    assert u.allclose(generic_map.heliographic_latitude, Latitude(sun.B0(generic_map.date)))


def test_heliographic_longitude(generic_map):
    # Needs a small tolerance to account for 32bit rounding errors
    assert u.allclose(generic_map.heliographic_longitude, 0 * u.deg, atol=1e-15*u.deg)


def test_units(generic_map):
    assert generic_map.spatial_units == ('arcsec', 'arcsec')


def test_cmap(generic_map):
    assert generic_map.cmap == matplotlib.colormaps['gray']


def test_coordinate_frame(aia171_test_map):
    frame = aia171_test_map.coordinate_frame
    assert isinstance(frame, sunpy.coordinates.Helioprojective)
    assert frame.observer.lat == aia171_test_map.observer_coordinate.frame.lat
    assert frame.observer.lon == aia171_test_map.observer_coordinate.frame.lon
    assert frame.observer.radius == aia171_test_map.observer_coordinate.frame.radius
    assert frame.obstime == aia171_test_map.reference_date


def test_heliographic_longitude_crln(hmi_test_map):
    assert_quantity_allclose(hmi_test_map.heliographic_longitude,
                             hmi_test_map.carrington_longitude - sun.L0(hmi_test_map.reference_date),
                             rtol=1e-3)  # A tolerance is needed because L0 is for Earth, not SDO


def test_observer_hgln_crln_priority():
    """
    When extracting the observer information from a FITS header, ensure
    Stonyhurst (HG) coordinates are preferred over Carrington (CR) if present
    (if not overridden by an instrument-specific `Map` subclass). Note that
    `Map` creates a custom FITS header with a sanitized observer location, so
    we also test a directly-instantiated `WCS` object in the coordinates
    module.
    """
    data = np.ones([6, 6], dtype=np.float64)
    header = {'CRVAL1': 0,
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
              'CTYPE1': 'HPLN-TAN',
              'CTYPE2': 'HPLT-TAN',
              'date-obs': '1970-01-01T00:00:00',
              'mjd-obs': 40587,
              'hglt_obs': 0,
              'hgln_obs': 0,
              'crlt_obs': 2,
              'crln_obs': 2,
              'dsun_obs': 10,
              'rsun_ref': 690000000}
    generic_map = sunpy.map.Map((data, header))

    c = generic_map.pixel_to_world(0*u.pix, 0*u.pix)
    assert c.observer.lon == 0 * u.deg
    # Note: don't test whether crlt or hglt is used---according to
    # coordinates.wcs_utils._set_wcs_aux_obs_coord, those are expected to
    # always be the same and so the same one is always used

    c = generic_map.wcs.pixel_to_world(0, 0)
    assert c.observer.lon == 0 * u.deg


def test_remove_observers(aia171_test_map):
    aia171_test_map._remove_existing_observer_location()
    with pytest.warns(SunpyMetadataWarning,
                      match='Missing metadata for observer: assuming Earth-based observer.*'):
        aia171_test_map.observer_coordinate


def test_partially_missing_observers(generic_map):
    generic_map.meta['hglt_obs'] = 0
    generic_map.meta['hgln_obs'] = 0
    generic_map.meta['crlt_obs'] = 0
    generic_map.meta['crln_obs'] = 0
    generic_map.meta.pop('dsun_obs')
    with pytest.warns(SunpyMetadataWarning,
                      match="Missing metadata for observer: assuming Earth-based observer.\n"
                            "For frame 'heliographic_stonyhurst' the following metadata is missing: dsun_obs\n"
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


_PC_KEYWORDS = ['PC1_1', 'PC1_2', 'PC2_1', 'PC2_2']
_CD_KEYWORDS = ['CD1_1', 'CD1_2', 'CD2_1', 'CD2_2']


@pytest.mark.parametrize('key', ['PC', 'CD'])
@pytest.mark.parametrize('i', [1, 2])
@pytest.mark.parametrize('j', [1, 2])
def test_rotation_matrix_defaults(generic_map, i, j, key):
    # Check that missing rotation keywords are set to correct defaults
    #
    # Relevant bit of the FITS standard:
    #
    # > PCi_j – [floating point; defaults: 1.0 when i = j, 0.0 otherwise]
    # > if any CDi_j keywords are present in the HDU, all other unspecified CDi_j keywords shall default to zero
    for keyword in _PC_KEYWORDS + _CD_KEYWORDS:
        if keyword in generic_map.meta:
            del generic_map.meta[keyword]

    keyword = f'{key}{i}_{j}'
    # Arbitrary number
    generic_map.meta[keyword] = 1.2
    if key == 'CD':
        expected = np.zeros((2, 2))
        expected[i-1, j-1] = 0.12
    elif key == 'PC':
        expected = np.eye(2)
        expected[i-1, j-1] = 1.2

    rot_mat = generic_map.rotation_matrix
    np.testing.assert_equal(rot_mat, expected)


def test_rotation_matrix_cd_cdelt():
    data = np.ones([6, 6], dtype=np.float64)
    header = {
        'CRVAL1': 0,
        'CRVAL2': 0,
        'CRPIX1': 5,
        'CRPIX2': 5,
        'CDELT1': 2,
        'CDELT2': 3,
        'CD1_1': 1,
        'CD1_2': -2,
        'CD2_1': 3,
        'CD2_2': 6,
        'NAXIS1': 6,
        'NAXIS2': 6,
        'CUNIT1': 'arcsec',
        'CUNIT2': 'arcsec',
        'CTYPE1': 'HPLN-TAN',
        'CTYPE2': 'HPLT-TAN',
    }
    cd_map = sunpy.map.Map((data, header))
    np.testing.assert_allclose(cd_map.rotation_matrix, np.array([[0.5, -1.], [1., 2.]]))


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
        'CUNIT2': 'arcsec',
        'CTYPE1': 'HPLN-TAN',
        'CTYPE2': 'HPLT-TAN',
    }
    cd_map = sunpy.map.Map((data, header))
    np.testing.assert_allclose(cd_map.rotation_matrix, np.array([[0., -1], [1., 0]]))


def test_swap_cd():
    amap = get_dummy_map_from_header(get_test_filepath('swap_lv1_20140606_000113.header'))
    np.testing.assert_allclose(amap.rotation_matrix, np.array([[1., 0], [0, 1.]]))


@pytest.mark.filterwarnings('ignore:Missing metadata for observer')
def test_crota_scale():
    # Test non-zero crota and unequal CDELT{1,2}
    n = 6
    data = np.ones([n, n], dtype=np.float64)
    header = {
        'CRVAL1': 0,
        'CRVAL2': 0,
        'CRPIX1': (n + 1) / 2,
        'CRPIX2': (n + 1) / 2,
        'NAXIS1': n,
        'NAXIS2': n,
        'CUNIT1': 'arcsec',
        'CUNIT2': 'arcsec',
        'CTYPE1': 'HPLN-TAN',
        'CTYPE2': 'HPLT-TAN',
        'DATE-OBS': '2020-01-01 00:00:00'
    }

    header.update({'CROTA2': 0, 'CDELT1': 1, 'CDELT2': 2})
    map1 = sunpy.map.Map(data, header)
    header.update({'CROTA2': 90, 'CDELT1': 2, 'CDELT2': 1})
    map2 = sunpy.map.Map(data, header)

    # Lower left coord
    coord1 = map1.pixel_to_world(*(-0.5, -0.5) * u.pix)
    # After rotating by 90 deg about the center of the map (CRPIX),
    # this should map to the lower right coordinate
    coord2 = map2.pixel_to_world(*(-0.5, n - 0.5) * u.pix)
    assert coord1.separation(coord2) < 1e-6 * u.arcsec


def test_world_to_pixel(generic_map):
    """Make sure conversion from data units to pixels is internally consistent"""
    test_pixel = generic_map.world_to_pixel(generic_map.reference_coordinate)
    assert_quantity_allclose(test_pixel, generic_map.reference_pixel)


def test_world_to_pixel_error(generic_map):
    strerr = 'Expected the following order of world arguments: SkyCoord'
    with pytest.raises(ValueError, match=strerr):
        generic_map.world_to_pixel(1)


def test_world_pixel_roundtrip(simple_map):
    pix = 1 * u.pix, 1 * u.pix
    coord = simple_map.pixel_to_world(*pix)
    pix_roundtrip = simple_map.world_to_pixel(coord)

    assert u.allclose(pix_roundtrip.x, pix[0], atol=1e-10 * u.pix)
    assert u.allclose(pix_roundtrip.y, pix[1], atol=1e-10 * u.pix)


def test_swapped_ctypes(simple_map):
    # Check that CTYPES different from normal work fine
    simple_map.meta['ctype1'] = 'HPLT-TAN'   # Usually HPLN
    simple_map.meta['ctype2'] = 'HPLN-TAN'   # Usually HPLT
    assert u.allclose(simple_map.bottom_left_coord.Tx, -4 * u.arcsec)
    assert u.allclose(simple_map.bottom_left_coord.Ty, -8 * u.arcsec)
    assert u.allclose(simple_map.top_right_coord.Tx, 4 * u.arcsec)
    assert u.allclose(simple_map.top_right_coord.Ty, 8 * u.arcsec)

    # Put them back
    simple_map.meta['ctype1'] = 'HPLN-TAN'   # Usually HPLN
    simple_map.meta['ctype2'] = 'HPLT-TAN'   # Usually HPLT
    assert u.allclose(simple_map.bottom_left_coord.Tx, -8 * u.arcsec)
    assert u.allclose(simple_map.bottom_left_coord.Ty, -4 * u.arcsec)
    assert u.allclose(simple_map.top_right_coord.Tx, 8 * u.arcsec)
    assert u.allclose(simple_map.top_right_coord.Ty, 4 * u.arcsec)


def test_save(aia171_test_map):
    """Tests the map save function"""
    aiamap = aia171_test_map
    afilename = tempfile.NamedTemporaryFile(suffix='fits').name
    aiamap.save(afilename, filetype='fits', overwrite=True)
    loaded_save = sunpy.map.Map(afilename)
    assert isinstance(loaded_save, sunpy.map.sources.AIAMap)
    # Compare metadata without considering ordering of keys
    assert loaded_save.meta.keys() == aiamap.meta.keys()
    for k in aiamap.meta:
        assert loaded_save.meta[k] == aiamap.meta[k]
    assert_quantity_allclose(loaded_save.data, aiamap.data)


@asdf_entry_points
def test_save_asdf(tmpdir, aia171_test_map):
    outpath = tmpdir/ "save_asdf.asdf"
    aia171_test_map.save(outpath, filetype= "asdf")
    loaded_save_asdf = sunpy.map.Map(str(outpath))
    assert isinstance(loaded_save_asdf, sunpy.map.sources.AIAMap)
    # Compare metadata without considering ordering of keys
    assert dict(loaded_save_asdf.meta) == dict(aia171_test_map.meta)
    np.testing.assert_array_equal(loaded_save_asdf.data, aia171_test_map.data)


def test_save_compressed(aia171_test_map):
    """Tests the map save function"""
    aiamap = aia171_test_map
    afilename = tempfile.NamedTemporaryFile(suffix='fits').name
    aiamap.save(afilename, filetype='fits', hdu_type=fits.CompImageHDU, overwrite=True)
    loaded_save = sunpy.map.Map(afilename)
    # We expect that round tripping to CompImageHDU will change the header and
    # the data a little.
    assert isinstance(loaded_save, sunpy.map.sources.AIAMap)


def test_shift_applied(generic_map):
    """Test that adding a shift actually updates the reference coordinate"""
    original_reference_coord = (generic_map.reference_coordinate.Tx,
                                generic_map.reference_coordinate.Ty)
    x_shift = 5 * u.arcsec
    y_shift = 13 * u.arcsec
    shifted_map = generic_map.shift_reference_coord(x_shift, y_shift)
    assert shifted_map.reference_coordinate.Tx - x_shift == original_reference_coord[0]
    assert shifted_map.reference_coordinate.Ty - y_shift == original_reference_coord[1]
    crval1 = ((generic_map.meta.get('crval1') * generic_map.spatial_units[0] +
               x_shift).to(shifted_map.spatial_units[0])).value
    assert shifted_map.meta.get('crval1') == crval1
    crval2 = ((generic_map.meta.get('crval2') * generic_map.spatial_units[1] +
               y_shift).to(shifted_map.spatial_units[1])).value
    assert shifted_map.meta.get('crval2') == crval2


def test_set_shift(generic_map):
    """Test that previously applied shift is stored in the shifted_value property"""
    x_shift = 5 * u.arcsec
    y_shift = 13 * u.arcsec
    shifted_map = generic_map.shift_reference_coord(x_shift, y_shift)
    mod_crval1 = shifted_map.meta.modified_items['crval1']
    mod_crval2 = shifted_map.meta.modified_items['crval2']
    assert x_shift == (mod_crval1.current - mod_crval1.original) * shifted_map.spatial_units[0]
    assert y_shift == (mod_crval2.current - mod_crval2.original) * shifted_map.spatial_units[1]


def test_shift_history(generic_map):
    """Test the shifted_value is added to a non-zero previous shift"""
    x_shift1 = 5 * u.arcsec
    y_shift1 = 13 * u.arcsec
    shifted_map1 = generic_map.shift_reference_coord(x_shift1, y_shift1)

    x_shift2 = -28.5 * u.arcsec
    y_shift2 = 120 * u.arcsec
    final_shifted_map = shifted_map1.shift_reference_coord(x_shift2, y_shift2)

    mod_crval1 = final_shifted_map.meta.modified_items['crval1']
    mod_crval2 = final_shifted_map.meta.modified_items['crval2']
    delta_crval1 = (mod_crval1.current - mod_crval1.original) * final_shifted_map.spatial_units[0]
    delta_crval2 = (mod_crval2.current - mod_crval2.original) * final_shifted_map.spatial_units[1]
    assert x_shift1 + x_shift2 == delta_crval1
    assert y_shift1 + y_shift2 == delta_crval2


def test_corners(simple_map):
    # These are the centers of the corner pixels
    assert u.allclose(simple_map.top_right_coord.Tx, 8 * u.arcsec)
    assert u.allclose(simple_map.top_right_coord.Ty, 4 * u.arcsec)
    assert u.allclose(simple_map.bottom_left_coord.Tx, -8 * u.arcsec)
    assert u.allclose(simple_map.bottom_left_coord.Ty, -4 * u.arcsec)


def test_center(simple_map):
    assert u.allclose(simple_map.center.Tx, 0 * u.arcsec, atol=1e-26 * u.arcsec)
    assert u.allclose(simple_map.center.Ty, 0 * u.arcsec)


def test_dimensions(simple_map):
    assert simple_map.dimensions[0] == 9 * u.pix
    assert simple_map.dimensions[1] == 9 * u.pix


pixel_corners = [
    [([0, 0] * u.pix, [0, 0] * u.pix), np.array([[0]])],
    [([-1, -1] * u.pix, [0, 0] * u.pix), np.array([[0]])],
    # 0.5, 0.5 is the edge of the first pixel, so make sure
    # we don't include any other pixels
    [([0, 0] * u.pix, [0.5, 0.5] * u.pix), np.array([[0]])],
    [([0, 0] * u.pix, [0, 0.51] * u.pix), np.array([[0],
                                                    [9]])],
    [([0, 0] * u.pix, [0.51, 0] * u.pix), np.array([[0, 1]])],
    [([0, 0] * u.pix, [0.51, 0.51] * u.pix), np.array([[0, 1],
                                                       [9, 10]])],
    [([0.1, 0.1] * u.pix, [1.6, 1.4] * u.pix), np.array([[0, 1, 2],
                                                         [9, 10, 11]])],
    [([0, 0] * u.pix, [20, 20] * u.pix), np.arange(81).reshape((9, 9))],
]
@pytest.mark.parametrize(("rect", "submap_out"), pixel_corners)
def test_submap_pixel(simple_map, rect, submap_out):
    # Check that result is the same specifying corners either way round
    for r in [dict(bottom_left=rect[0], top_right=rect[1]),
              dict(bottom_left=rect[1], top_right=rect[0])]:
        submap = simple_map.submap(**r)
        np.testing.assert_equal(submap.data, submap_out)


# The (0.5, 0.5) case is skipped as boundary points cannot reliably tested when
# converting to world coordinates due to round-off error when round-tripping
# through pixel_to_world -> world_to_pixel
@pytest.mark.parametrize(("rect", "submap_out"), pixel_corners[:2] + pixel_corners[3:])
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


@pytest.mark.parametrize('test_map', ["aia171_roll_map", "aia171_test_map",
                                      "hmi_test_map", "aia171_test_map_with_mask"],
                         indirect=['test_map'])
def test_submap_world_corners(test_map):
    """
    This test checks that when an unaligned map is cropped with submap that the
    resulting map contains all four corners of the input world coordinate
    bounding box.
    """
    corners = SkyCoord(Tx=[300, 300, 800, 800], Ty=[0, 500, 500, 0],
                       unit=u.arcsec, frame=test_map.coordinate_frame)

    submap = test_map.submap(corners[0], top_right=corners[2])

    pix_corners = np.array(submap.wcs.world_to_pixel(corners)).T
    for pix_corner in pix_corners:
        assert ((-0.5, -0.5) <= pix_corner).all()
        assert (pix_corner <= submap.data.shape[::-1]).all()

    if test_map.mask is not None:
        assert submap.mask.shape == submap.data.shape


@pytest.mark.parametrize('test_map', ["aia171_test_map", "heliographic_test_map"],
                         indirect=['test_map'])
def test_submap_hgs_corners(test_map):
    """
    This test checks that when an unaligned map is cropped with submap that the
    resulting map contains all four corners of the input world coordinate
    bounding box.
    """
    corners = SkyCoord([10, 10, 40, 40], [-10, 30, 30, -10],
                       unit=u.deg, frame="heliographic_stonyhurst",
                       obstime=test_map.date)

    submap = test_map.submap(corners[0], top_right=corners[2])

    pix_corners = np.array(submap.wcs.world_to_pixel(corners)).T
    for pix_corner in pix_corners:
        assert ((-0.5, -0.5) <= pix_corner).all()
        assert (pix_corner <= submap.data.shape[::-1]).all()


# Check that submap works with units convertible to pix but that aren't pix
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
    assert simple_map.reference_pixel.x == 4 * u.pix
    assert simple_map.reference_pixel.y == 4 * u.pix


@pytest.mark.parametrize('shape', [[1, 1], [3, 3]])
def test_resample(simple_map, shape):
    resampled = simple_map.resample(shape * u.pix, method='linear')
    assert np.mean(resampled.data) == np.mean(simple_map.data)
    if shape == [1, 1]:
        # Should be the mean of [0,1,2,...78,79,80]
        assert resampled.data == np.array([[40]])

    # Check that the corner coordinates of the input and output are the same
    resampled_lower_left = resampled.pixel_to_world(-0.5 * u.pix, -0.5 * u.pix)
    original_lower_left = simple_map.pixel_to_world(-0.5 * u.pix, -0.5 * u.pix)
    assert u.allclose(resampled_lower_left.Tx, original_lower_left.Tx)
    assert u.allclose(resampled_lower_left.Ty, original_lower_left.Ty)

    resampled_upper_right = resampled.pixel_to_world((shape[0] - 0.5) * u.pix,
                                                    (shape[1] - 0.5) * u.pix)
    original_upper_right = simple_map.pixel_to_world(8.5 * u.pix, 8.5 * u.pix)
    assert u.allclose(resampled_upper_right.Tx, original_upper_right.Tx)
    assert u.allclose(resampled_upper_right.Ty, original_upper_right.Ty)


resample_test_data = [('linear', (100, 200) * u.pixel),
                      ('nearest', (512, 128) * u.pixel),
                      ('spline', (200, 200) * u.pixel)]


@pytest.mark.parametrize(("sample_method", "new_dimensions"), resample_test_data)
def test_resample_dimensions(generic_map, sample_method, new_dimensions):
    """Check that resampled map has expected dimensions."""
    resampled_map = generic_map.resample(new_dimensions, method=sample_method)
    assert resampled_map.dimensions[0] == new_dimensions[0]
    assert resampled_map.dimensions[1] == new_dimensions[1]


@pytest.mark.parametrize(("sample_method", "new_dimensions"), resample_test_data)
def test_resample_metadata(generic_map, sample_method, new_dimensions):
    """
    Check that the resampled map has correctly adjusted metadata.
    """
    resampled_map = generic_map.resample(new_dimensions, method=sample_method)

    scale_x = generic_map.data.shape[1] / resampled_map.data.shape[1]
    scale_y = generic_map.data.shape[0] / resampled_map.data.shape[0]

    assert resampled_map.meta['cdelt1'] == scale_x * generic_map.meta['cdelt1']
    assert resampled_map.meta['cdelt2'] == scale_y * generic_map.meta['cdelt2']
    assert resampled_map.meta['pc1_1'] == generic_map.meta['pc1_1']
    assert resampled_map.meta['pc1_2'] == scale_y / scale_x * generic_map.meta['pc1_2']
    assert resampled_map.meta['pc2_1'] == scale_x / scale_y * generic_map.meta['pc2_1']
    assert resampled_map.meta['pc2_2'] == generic_map.meta['pc2_2']

    # TODO: we should really test the numbers here, not just that the correct
    # header values have been modified. However, I am lazy and we have figure
    # tests.
    assert resampled_map.meta['crpix1'] != generic_map.meta['crpix1']
    assert resampled_map.meta['crpix2'] != generic_map.meta['crpix2']
    assert u.allclose(resampled_map.meta['crval1'], generic_map.meta['crval1'])
    assert u.allclose(resampled_map.meta['crval2'], generic_map.meta['crval2'])
    assert resampled_map.meta['naxis1'] == new_dimensions[0].value
    assert resampled_map.meta['naxis2'] == new_dimensions[1].value
    for key in generic_map.meta:
        if key not in ('cdelt1', 'cdelt2', 'pc1_2', 'pc2_1', 'crpix1', 'crpix2', 'crval1',
                       'crval2', 'naxis1', 'naxis2'):
            assert resampled_map.meta[key] == generic_map.meta[key]


@pytest.mark.parametrize(("sample_method", "new_dimensions"), resample_test_data)
def test_resample_simple_map(simple_map, sample_method, new_dimensions):
    # Put the reference pixel at the top-right of the bottom-left pixel
    simple_map.meta['crpix1'] = 1.5
    simple_map.meta['crpix2'] = 1.5
    assert list(simple_map.reference_pixel) == [0.5 * u.pix, 0.5 * u.pix]
    # Make the superpixel map
    new_dims = (9, 6) * u.pix
    resamp_map = simple_map.resample(new_dims, method=sample_method)
    # Reference pixel should change, but reference coordinate should not
    assert u.allclose(list(resamp_map.reference_pixel), [0.5 * u.pix, 0.16666667 * u.pix])
    assert resamp_map.reference_coordinate == simple_map.reference_coordinate


def test_superpixel_simple_map(simple_map):
    # Put the reference pixel at the top-right of the bottom-left pixel
    simple_map.meta['crpix1'] = 1.5
    simple_map.meta['crpix2'] = 1.5
    assert list(simple_map.reference_pixel) == [0.5 * u.pix, 0.5 * u.pix]
    # Make the superpixel map
    new_dims = (2, 2) * u.pix
    superpix_map = simple_map.superpixel(new_dims)
    # Reference pixel should change, but reference coordinate should not
    assert list(superpix_map.reference_pixel) == [0 * u.pix, 0 * u.pix]
    assert superpix_map.reference_coordinate == simple_map.reference_coordinate

    # Check that offset works
    superpix_map = simple_map.superpixel(new_dims, offset=[1, 2] * u.pix)
    # Reference pixel should change, but reference coordinate should not
    assert u.allclose(list(superpix_map.reference_pixel),
                      [-0.5 * u.pix, -1 * u.pix])
    assert superpix_map.reference_coordinate == simple_map.reference_coordinate


@pytest.mark.parametrize('f', [np.sum, np.mean])
def test_superpixel_dims_values(aia171_test_map, f):
    dimensions = (2, 2) * u.pix
    superpix_map = aia171_test_map.superpixel(dimensions, func=f)

    # Check dimensions of new map
    old_dims = u.Quantity(aia171_test_map.dimensions)
    expected_new_dims = old_dims * (1 * u.pix / dimensions)
    assert superpix_map.dimensions[0] == expected_new_dims[0]
    assert superpix_map.dimensions[1] == expected_new_dims[1]

    # Check value of lower left pixel is calculated correctly
    expected = f(aia171_test_map.data[0:2, 0:2])
    assert_quantity_allclose(superpix_map.data[0, 0], expected)


@pytest.mark.parametrize(("f", "dimensions"), [(np.sum, (2, 3)*u.pix),
                                               (np.mean, (3, 2)*u.pix)])
def test_superpixel_metadata(generic_map, f, dimensions):
    superpix_map = generic_map.superpixel(dimensions, func=f)

    scale_x, scale_y = dimensions.value

    assert superpix_map.meta['cdelt1'] == scale_x * generic_map.meta['cdelt1']
    assert superpix_map.meta['cdelt2'] == scale_y * generic_map.meta['cdelt2']
    assert superpix_map.meta['pc1_1'] == generic_map.meta['pc1_1']
    assert superpix_map.meta['pc1_2'] == scale_y / scale_x * generic_map.meta['pc1_2']
    assert superpix_map.meta['pc2_1'] == scale_x / scale_y * generic_map.meta['pc2_1']
    assert superpix_map.meta['pc2_2'] == generic_map.meta['pc2_2']

    assert superpix_map.meta['crpix1'] - 0.5 == (generic_map.meta['crpix1'] - 0.5) / scale_x
    assert superpix_map.meta['crpix2'] - 0.5 == (generic_map.meta['crpix2'] - 0.5) / scale_y
    assert u.allclose(superpix_map.meta['crval1'], generic_map.meta['crval1'])
    assert u.allclose(superpix_map.meta['crval2'], generic_map.meta['crval2'])
    assert superpix_map.meta['naxis1'] == generic_map.meta['naxis1'] / scale_x
    assert superpix_map.meta['naxis2'] == generic_map.meta['naxis2'] / scale_y
    for key in generic_map.meta:
        if key not in ('cdelt1', 'cdelt2', 'pc1_2', 'pc2_1', 'crpix1', 'crpix2', 'crval1',
                       'crval2', 'naxis1', 'naxis2'):
            assert superpix_map.meta[key] == generic_map.meta[key]


def test_superpixel_masked(aia171_test_map_with_mask):
    input_dims = u.Quantity(aia171_test_map_with_mask.dimensions)
    dimensions = (2, 2) * u.pix
    # Test that the mask is respected
    superpix_map = aia171_test_map_with_mask.superpixel(dimensions)
    assert superpix_map.mask is not None
    # Check the shape of the mask
    expected_shape = input_dims * (1 * u.pix / dimensions)
    assert np.all(superpix_map.mask.shape * u.pix == expected_shape)

    # Test that the offset is respected
    superpix_map = aia171_test_map_with_mask.superpixel(dimensions, offset=(1, 1) * u.pix)
    assert superpix_map.dimensions[0] == expected_shape[0] - 1 * u.pix
    assert superpix_map.dimensions[1] == expected_shape[1] - 1 * u.pix

    dimensions = (7, 9) * u.pix
    superpix_map = aia171_test_map_with_mask.superpixel(dimensions, offset=(4, 4) * u.pix)
    expected_shape = np.round(input_dims * (1 * u.pix / dimensions))
    assert superpix_map.dimensions[0] == expected_shape[0] - 1 * u.pix
    assert superpix_map.dimensions[1] == expected_shape[1] - 1 * u.pix


def test_superpixel_masked_conservative_mask_true(aia171_test_map_with_mask):
    input_dims = u.Quantity(aia171_test_map_with_mask.dimensions)
    dimensions = (2, 2) * u.pix

    superpix_map = aia171_test_map_with_mask.superpixel(dimensions, conservative_mask=True)
    assert superpix_map.mask is not None

    expected_shape = input_dims * (1 * u.pix / dimensions)
    assert np.all(superpix_map.mask.shape * u.pix == expected_shape)

    # Verify mask values (bin_mask=True)
    reshaped_mask = reshape_image_to_4d_superpixel(
        aia171_test_map_with_mask.mask,
        [dimensions[1].value, dimensions[0].value],
        [0, 0],
    )
    expected_mask = np.any(reshaped_mask, axis=(1, 3))
    assert np.array_equal(superpix_map.mask, expected_mask)


def test_superpixel_units(generic_map):
    new_dims = (2, 2) * u.pix
    super1 = generic_map.superpixel(new_dims)
    super2 = generic_map.superpixel(new_dims.to(u.kpix))
    assert super1.meta == super2.meta

    offset = (1, 2) * u.pix
    super1 = generic_map.superpixel(new_dims, offset=offset)
    super2 = generic_map.superpixel(new_dims, offset=offset.to(u.kpix))
    assert super1.meta == super2.meta


def test_superpixel_fractional_inputs(generic_map):
    super1 = generic_map.superpixel((2, 3) * u.pix)
    super2 = generic_map.superpixel((2.2, 3.2) * u.pix)
    assert np.all(super1.data == super2.data)
    assert super1.meta == super2.meta


@pytest.mark.parametrize('method', ['resample', 'superpixel'])
@settings(
    max_examples=10,
    # Lots of draws can be discarded when checking matrix is non-singular
    suppress_health_check=[HealthCheck.filter_too_much, HealthCheck.function_scoped_fixture],
    deadline=1000,
)
@given(pc=matrix_meta('pc'))
def test_resample_rotated_map_pc(pc, method, simple_map):
    smap = deepcopy(simple_map)
    smap.meta.update(pc)
    # Check superpixel with a rotated map with unequal resampling
    new_dims = (1, 2) * u.pix
    new_map = getattr(smap, method)(new_dims)
    # Coordinate of the lower left corner should not change
    ll_pix = [-0.5, -0.5]*u.pix
    assert smap.pixel_to_world(*ll_pix).separation(
        new_map.pixel_to_world(*ll_pix)).to(u.arcsec) < 1e-8 * u.arcsec


@pytest.mark.parametrize('method', ['resample', 'superpixel'])
@settings(
    max_examples=10,
    # Lots of draws can be discarded when checking matrix is non-singular
    suppress_health_check=[HealthCheck.filter_too_much, HealthCheck.function_scoped_fixture],
    deadline=1000,
)
@given(cd=matrix_meta('cd'))
def test_resample_rotated_map_cd(cd, method, simple_map):
    smap = deepcopy(simple_map)
    smap.meta.update(cd)
    for key in ['cdelt1', 'cdelt2', 'pc1_1', 'pc1_2', 'pc2_1', 'pc2_2']:
        del smap.meta[key]
    # Check superpixel with a rotated map with unequal resampling
    new_dims = (1, 2) * u.pix
    new_map = getattr(smap, method)(new_dims)
    # Coordinate of the lower left corner should not change
    ll_pix = [-0.5, -0.5]*u.pix
    assert smap.pixel_to_world(*ll_pix).separation(
        new_map.pixel_to_world(*ll_pix)).to(u.arcsec) < 1e-8 * u.arcsec


def test_superpixel_err(generic_map):
    with pytest.raises(ValueError, match="Offset is strictly non-negative."):
        generic_map.superpixel((2, 2) * u.pix, offset=(-2, 2) * u.pix)


def calc_new_matrix(angle):
    c = np.cos(np.deg2rad(angle))
    s = np.sin(np.deg2rad(angle))
    return np.array([[c, -s], [s, c]])


def test_rotate(aia171_test_map):
    # We use order=0 for many of these tests to minimize losing edge pixels due to interpolation
    # with NaNs that are used as the default `missing` value

    rotated_map_1 = aia171_test_map.rotate(20 * u.deg, order=0)
    rotated_map_2 = rotated_map_1.rotate(20 * u.deg, order=0)
    np.testing.assert_allclose(rotated_map_1.rotation_matrix,
                               np.dot(aia171_test_map.rotation_matrix, calc_new_matrix(20).T))
    np.testing.assert_allclose(rotated_map_2.rotation_matrix,
                               np.dot(aia171_test_map.rotation_matrix, calc_new_matrix(40).T))

    # Rotation of a map by a non-integral multiple of 90 degrees expands the map
    # and assigns the value of NaN to corner regions. The mean will be approximately
    # the same, although there will be slight change due to the loss of edge pixels
    # due to interpolation with the NaNs.
    assert rotated_map_2.data.shape > rotated_map_1.data.shape > aia171_test_map.data.shape
    assert np.isnan(rotated_map_1.data[0, 0])
    assert np.isnan(rotated_map_2.data[0, 0])
    np.testing.assert_allclose(aia171_test_map.mean(), rotated_map_1.mean(), rtol=5e-3)
    np.testing.assert_allclose(aia171_test_map.mean(), rotated_map_2.mean(), rtol=5e-3)

    # A scaled-up map should have the same mean because the output map should be expanded
    rotated_map_3 = aia171_test_map.rotate(0 * u.deg, order=0, scale=2)
    np.testing.assert_allclose(aia171_test_map.mean(), rotated_map_3.mean(), rtol=1e-4)

    # Mean and std should be equal for a 90 degree rotation as long as 1 pixel is cropped out on
    # all sides
    rotated_map_4 = aia171_test_map.rotate(90 * u.deg, order=0)
    np.testing.assert_allclose(aia171_test_map.data[1:-1, 1:-1].mean(),
                               rotated_map_4.data[1:-1, 1:-1].mean(), rtol=1e-10)
    np.testing.assert_allclose(aia171_test_map.data[1:-1, 1:-1].std(),
                               rotated_map_4.data[1:-1, 1:-1].std(), rtol=1e-10)

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


def test_rotate_with_incompatible_missing_dtype_error():
    data = np.arange(0, 100).reshape(10, 10)
    coord = SkyCoord(0 * u.arcsec, 0 * u.arcsec, obstime='2013-10-28',
                     observer='earth', frame=sunpy.coordinates.Helioprojective)
    header = sunpy.map.make_fitswcs_header(data, coord)
    test_map = sunpy.map.Map(data, header)
    with pytest.raises(ValueError, match="The underlying data is integers, but the fill value for "
                                         "missing pixels cannot be cast to an integer"):
        test_map.rotate(angle=45 * u.deg, missing=np.nan, order=3)


def test_rotate_crpix_zero_degrees(generic_map):
    # Rotating a map by zero degrees should not change the location of the reference pixel at all
    rotated_map = generic_map.rotate(0*u.deg)
    assert rotated_map.reference_pixel.x == generic_map.reference_pixel.x
    assert rotated_map.reference_pixel.y == generic_map.reference_pixel.y


def test_rotate_pad_crpix(generic_map):
    rotated_map = generic_map.rotate(30*u.deg)
    # This tests that the reference pixel of the map is in the expected place.
    assert rotated_map.data.shape != generic_map.data.shape
    assert_quantity_allclose(u.Quantity(rotated_map.reference_pixel),
                             u.Quantity((5.04903811, 6.54903811), u.pix))


def test_rotate_recenter(generic_map):
    rotated_map = generic_map.rotate(20 * u.deg, recenter=True)
    pixel_array_center = (np.flipud(rotated_map.data.shape) - 1) / 2.0

    assert_quantity_allclose(
        pixel_array_center * u.pix, u.Quantity(rotated_map.reference_pixel))


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
    rot_map = generic_map.rotate(rmatrix=np.array([[0, 1], [-1, 0]]))
    np.testing.assert_allclose(rot_map.rotation_matrix, np.array([[-1, 0], [0, -1]]))


def test_rotate_rmatrix_angle(generic_map):
    with pytest.raises(ValueError, match="You cannot specify both an angle and a rotation matrix."):
        generic_map.rotate(angle=5*u.deg, rmatrix=np.array([[1, 0], [0, 1]]))


def test_rotate_invalid_order(generic_map):
    with pytest.raises(ValueError, match="Order must be between 0 and 5."):
        generic_map.rotate(order=6)
    with pytest.raises(ValueError, match="Order must be between 0 and 5."):
        generic_map.rotate(order=-1)


def test_rotate_assumed_obstime():
    # Create an HPC map that is missing the observing time and has an off-disk reference coordinate
    header = {
        'crval1': -2000,
        'crval2': 0,
        'cdelt1': 1,
        'cdelt2': 1,
        'ctype1': 'HPLN-TAN',
        'ctype2': 'HPLT-TAN',
        'naxis': 2,
        'naxis1': 10,
        'naxis2': 10,
        'cunit1': 'arcsec',
        'cunit2': 'arcsec',
        'crpix1': 4.5,
        'crpix2': 6.5,
        'hglt_obs': 0,
        'hgln_obs': 0,
        'dsun_obs': 150000000000,
        'rsun_ref': 700000000,
    }
    original = sunpy.map.Map(np.zeros((10, 10)), header)

    # Accessing the date makes the assumption of "now" for obstime
    with pytest.warns(SunpyMetadataWarning, match="Missing metadata for observation time"):
        original.date

    # The assumption has already been made, so no further warning should be emitted by rotate()
    rotated = original.rotate(0*u.deg)

    # The reference coordinate should be unchanged by this 0-degree rotation
    # Since the reference coordinate is off-disk, a non-identity transformation would result in NaNs
    assert_quantity_allclose(rotated.reference_pixel.x, original.reference_pixel.x)
    assert_quantity_allclose(rotated.reference_pixel.y, original.reference_pixel.y)

    # The returned map should also be missing observing time
    with pytest.warns(SunpyMetadataWarning, match="Missing metadata for observation time"):
        rotated.date


def test_as_mpl_axes_aia171(aia171_test_map):
    ax = plt.subplot(projection=aia171_test_map)
    assert isinstance(ax, wcsaxes.WCSAxes)
    assert all([ct1 == ct2 for ct1, ct2 in zip(ax.wcs.wcs.ctype, aia171_test_map.wcs.wcs.ctype)])


def test_plot_with_norm_none(aia171_test_map):
    # Confirm that norm == None does not raise an error, see https://github.com/sunpy/sunpy/pull/7261
    ax = plt.subplot(projection=aia171_test_map)
    aia171_test_map.plot(axes=ax, norm=None, vmin=0, vmax=0)


def test_validate_meta(generic_map):
    """Check to see if_validate_meta displays an appropriate error"""
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
    with pytest.warns(SunpyMetadataWarning) as w:
        sunpy.map.Map((generic_map.data, bad_header))

    assert 'waveunit'.upper() in str(w[0].message)


def test_validate_non_spatial(generic_map):
    generic_map.meta['cunit2'] = 'Angstrom'
    err_msg = ("Map only supports spherical coordinate systems with angular units "
               "(ie. equivalent to arcsec), but this map has units ['arcsec', 'Angstrom']")
    with pytest.raises(sunpy.map.MapMetaValidationError, match=re.escape(err_msg)):
        sunpy.map.Map(generic_map.data, generic_map.meta)


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


@pytest.mark.skipif(pytest.__version__ < "8.0.0", reason="pytest >= 8.0.0 raises two warnings for this test")
def test_more_than_two_dimensions():
    """
    Checks to see if an appropriate error is raised when a FITS with more than two dimensions is
    loaded. We need to load a >2-dim dataset with a TELESCOP header
    """

    # Data crudely represents 4 stokes, 4 wavelengths with Y,X of 3 and 5.
    bad_data = np.random.rand(4, 4, 3, 5)
    hdr = fits.Header()
    hdr['TELESCOP'] = 'XXX'
    hdr['cunit1'] = 'arcsec'
    hdr['cunit2'] = 'arcsec'
    with pytest.warns(SunpyMetadataWarning, match='Missing CTYPE'):
        with pytest.warns(SunpyUserWarning, match='This file contains more than 2 dimensions.'):
            bad_map = sunpy.map.Map(bad_data, hdr)
    # Test fails if map.ndim > 2 and if the dimensions of the array are wrong.
    assert bad_map.ndim == 2
    assert_quantity_allclose(bad_map.dimensions, (5, 3) * u.pix)


def test_missing_metadata_warnings():
    # Checks that warnings for missing metadata are only raised once
    header = {
        'cunit1': 'arcsec',
        'cunit2': 'arcsec',
        'ctype1': 'HPLN-TAN',
        'ctype2': 'HPLT-TAN',
    }
    with pytest.warns(Warning) as record:  # NOQA: PT030,PT031
        array_map = sunpy.map.Map(np.random.rand(20, 15), header)
        array_map.peek()
    # There should be 2 warnings for missing metadata (obstime and observer location)
    assert len([w for w in record if w.category in (SunpyMetadataWarning, SunpyUserWarning)]) == 2


def test_fits_header(aia171_test_map):
    assert isinstance(aia171_test_map.fits_header, fits.Header)


def test_bad_coordframe_repr(generic_map):
    generic_map.meta['CTYPE1'] = "STUART1"
    generic_map.meta['CTYPE2'] = "STUART2"
    with pytest.warns(UserWarning,
                      match="Could not determine coordinate frame from map metadata"):
        assert 'Unknown' in generic_map.__repr__()


def test_non_str_key():
    header = {'cunit1': 'arcsec',
              'cunit2': 'arcsec',
              None: None,  # Cannot parse this into WCS
              }
    with pytest.raises(ValueError, match='All MetaDict keys must be strings'):
        sunpy.map.GenericMap(np.zeros((10, 10)), header)


def test_updating_of_naxisi_on_rotate(aia171_test_map):
    aia171_test_map_rotated = aia171_test_map.rotate(45 * u.deg, missing=0)
    assert aia171_test_map.data.shape == (128, 128)
    assert aia171_test_map_rotated.meta['NAXIS1'] == 182
    assert aia171_test_map_rotated.meta['NAXIS2'] == 182


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


def test_quicklook(mocker, aia171_test_map):
    mockwbopen = mocker.patch('webbrowser.open_new_tab')
    aia171_test_map.quicklook()
    # Check that the mock web browser was opened with a file URL
    mockwbopen.assert_called_once()
    file_url = mockwbopen.call_args[0][0]
    assert file_url.startswith('file://')
    # Open the file specified in the URL and confirm that it contains the HTML
    with open(file_url[7:]) as f:
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


def test_submap_inputs(generic_map2, coords):
    bl_coord, tr_coord, bl_tr_coord = coords

    inputs = (
        ((bl_coord,), dict(top_right=tr_coord)),
        ((bl_coord,), dict(width=width_deg, height=height_deg)),
        ((bl_tr_coord,), {}),
        ((bl_pix,), dict(top_right=tr_pix)),
        ((bl_pix,), dict(width=width_pix, height=height_pix)),
        ((bl_tr_coord.frame,), {}),
    )

    for args, kwargs in inputs:
        smap = generic_map2.submap(*args, **kwargs)
        assert u.allclose(smap.dimensions, (3, 3) * u.pix)


def test_contour_deprecation_warning(simple_map):

    with pytest.warns(SunpyDeprecationWarning, match="The contour function is deprecated and may be removed in a future version.\\s+Use sunpy.map.GenericMap.find_contours instead."):
        simple_map.contour(1.5)


def test_find_contours_contourpy(simple_map):
    data = np.ones(simple_map.data.shape)
    data[4, 4] = 2
    simple_map = sunpy.map.Map(data, simple_map.meta)
    # 4 is the central pixel of the map, so contour half way between 1 and 2
    contours = simple_map.find_contours(1.5, method='contourpy')
    assert len(contours) == 1
    contour = contours[0]
    assert contour.observer.lat == simple_map.observer_coordinate.frame.lat
    assert contour.observer.lon == simple_map.observer_coordinate.frame.lon
    assert contour.obstime == simple_map.date
    assert u.allclose(contour.Tx, [-1, 0, 1, 0, -1] * u.arcsec, atol=1e-10 * u.arcsec)
    assert u.allclose(contour.Ty, [ 0, -0.5, 0, 0.5, 0] * u.arcsec, atol=1e-10 * u.arcsec)
    with pytest.raises(ValueError, match='level must be a single scalar value'):
        simple_map.find_contours([1.5, 2.5])


def test_find_contours_skimage(simple_map):
    data = np.ones(simple_map.data.shape)
    data[4, 4] = 2
    simple_map = sunpy.map.Map(data, simple_map.meta)
    # 4 is the central pixel of the map, so contour half way between 1 and 2
    contours = simple_map.find_contours(1.5, method='skimage')
    assert len(contours) == 1
    contour = contours[0]
    assert contour.observer.lat == simple_map.observer_coordinate.frame.lat
    assert contour.observer.lon == simple_map.observer_coordinate.frame.lon
    assert contour.obstime == simple_map.date
    assert u.allclose(contour.Tx, [0, -1, 0, 1, 0] * u.arcsec, atol=1e-10 * u.arcsec)
    assert u.allclose(contour.Ty, [0.5, 0, -0.5, 0, 0.5] * u.arcsec, atol=1e-10 * u.arcsec)
    with pytest.raises(ValueError, match='level must be a single scalar value'):
        simple_map.find_contours([1.5, 2.5])


def test_find_contours_invalid_library(simple_map):
    with pytest.raises(ValueError, match="Unknown method 'invalid_method'. Use 'contourpy' or 'skimage'."):
        simple_map.find_contours(1.5, method='invalid_method')


def test_find_contours_units(simple_map):
    # Check that contouring with units works as intended
    simple_map.meta['bunit'] = 'm'
    # Same units
    contours = simple_map.find_contours(1.5 * u.m)
    assert len(contours) == 1

    # Different units, but convertible
    contours_cm = simple_map.find_contours(150 * u.cm)
    for c1, c2 in zip(contours, contours_cm):
        assert np.all(c1 == c2)

    # Percentage
    contours_percent = simple_map.find_contours(50 * u.percent)
    high = np.max(simple_map.data)
    low = np.min(simple_map.data)
    middle = high - (high - low) / 2
    contours_ref = simple_map.find_contours(middle * simple_map.unit)
    for c1, c2 in zip(contours_percent, contours_ref):
        assert np.all(c1 == c2)


def test_find_contours_inputs(simple_map):
    with pytest.raises(ValueError, match='Contour levels must be increasing'):
        simple_map.draw_contours([10, -10] * u.dimensionless_unscaled)
    with pytest.raises(ValueError, match=re.escape('The provided level (1000.0) is not smaller than the maximum data value (80)')):
        simple_map.draw_contours(1000 * u.dimensionless_unscaled, fill=True)

    simple_map.meta['bunit'] = 'm'

    with pytest.raises(TypeError, match='The levels argument has no unit attribute'):
        simple_map.draw_contours(1.5)
    with pytest.raises(TypeError, match='The levels argument has no unit attribute'):
        simple_map.find_contours(1.5)

    with pytest.raises(u.UnitsError, match=re.escape("'s' (time) and 'm' (length) are not convertible")):
        simple_map.draw_contours(1.5 * u.s)
    with pytest.raises(u.UnitsError, match=re.escape("'s' (time) and 'm' (length) are not convertible")):
        simple_map.find_contours(1.5 * u.s)

    # With no units, check that dimensionless works
    simple_map.meta.pop('bunit')
    simple_map.draw_contours(1.5 * u.dimensionless_unscaled)
    simple_map.find_contours(1.5 * u.dimensionless_unscaled)

    with pytest.raises(u.UnitsError, match='This map has no unit'):
        simple_map.draw_contours(1.5 * u.m)
    with pytest.raises(u.UnitsError, match='This map has no unit'):
        simple_map.find_contours(1.5 * u.m)


def test_print_map(generic_map):
    out_repr = generic_map.__repr__()
    assert isinstance(out_repr, str)
    assert object.__repr__(generic_map) in out_repr
    out_str = generic_map.__str__()
    assert isinstance(out_str, str)
    assert out_str in out_repr


def test_parse_submap_quantity_inputs(aia171_test_map):
    bottom_left = (0, 0)*u.arcsec
    top_right = (200, 200)*u.arcsec
    width = 200*u.arcsec
    height = 300*u.arcsec

    with pytest.raises(ValueError, match=re.escape("Either top_right alone or both width and height "
                       "must be specified when bottom_left is a Quantity")):
        aia171_test_map.submap(bottom_left=bottom_left[0],
                               top_right=None, width=None, height=None)

    with pytest.raises(ValueError, match=re.escape("bottom_left must have shape (2, ) "
                       "when specified as a Quantity")):
        aia171_test_map.submap(bottom_left=bottom_left[0],
                               top_right=top_right, width=None, height=None)

    with pytest.raises(ValueError, match=re.escape("top_right must have shape (2, ) when specified as "
                       "a Quantity")):
        aia171_test_map.submap(bottom_left=bottom_left,
                               top_right=top_right[0], width=None, height=None)

    with pytest.raises(TypeError, match=re.escape("When bottom_left is a Quantity, top_right "
                       "must be a Quantity in units of pixels.")):
        aia171_test_map.submap(bottom_left=bottom_left,
                               top_right=top_right, width=None, height=None)

    with pytest.raises(TypeError, match=re.escape("When bottom_left is a Quantity, width and height "
                       "must be a Quantity in units of pixels.")):
        aia171_test_map.submap(bottom_left=bottom_left,
                               top_right=None, width=width, height=height)


def test_wavelength_properties(simple_map):
    simple_map.meta.pop('waveunit', None)
    simple_map.meta['wavelnth'] = 1
    assert simple_map.measurement == 1 * u.one
    assert simple_map.wavelength == 1 * u.one

    simple_map.meta['waveunit'] = ''
    assert simple_map.measurement == 1 * u.one
    assert simple_map.wavelength == 1 * u.one

    simple_map.meta['waveunit'] = 'm'
    assert simple_map.measurement == 1 * u.m
    assert simple_map.wavelength == 1 * u.m


def test_meta_modifications(aia171_test_map):
    aiamap = aia171_test_map
    old_cdelt1 = aiamap.meta['cdelt1']
    aiamap.meta['cdelt1'] = 20

    assert aiamap.meta.original_meta != aiamap.meta
    assert aiamap.meta.added_items == {}
    assert aiamap.meta.removed_items == {}
    assert aiamap.meta.modified_items == {'cdelt1': ModifiedItem(old_cdelt1, 20)}

    # Check that rotate doesn't modify the original metadata
    aiamap_rot = aiamap.rotate(30 * u.deg)
    assert aiamap_rot.meta.original_meta == aiamap.meta.original_meta
    assert set(aiamap_rot.meta.added_items.keys()) == set(['pc1_1', 'pc1_2', 'pc2_1', 'pc2_2'])
    assert set(aiamap_rot.meta.removed_items.keys()) == set(['crota2'])
    assert set(aiamap_rot.meta.modified_items) == set(['cdelt1', 'crpix1', 'crpix2', 'crval1', 'naxis1', 'naxis2'])


def test_no_wcs_observer_info(heliographic_test_map):
    # Check that HeliographicCarrington WCS has observer info set
    assert isinstance(heliographic_test_map.coordinate_frame, HeliographicCarrington)
    wcs_aux = heliographic_test_map.wcs.wcs.aux
    assert wcs_aux.hgln_obs is not None
    assert wcs_aux.hglt_obs is not None
    assert wcs_aux.dsun_obs is not None

    # Remove observer information, and change coordinate system to HeliographicStonyhurst
    heliographic_test_map.meta.pop('HGLN_OBS')
    heliographic_test_map.meta.pop('HGLT_OBS')
    heliographic_test_map.meta.pop('DSUN_OBS')
    heliographic_test_map.meta['CTYPE1'] = 'HGLN-CAR'
    heliographic_test_map.meta['CTYPE2'] = 'HGLT-CAR'
    assert isinstance(heliographic_test_map.coordinate_frame, HeliographicStonyhurst)

    # Check that GenericMap.wcs doesn't set an observer
    wcs_aux = heliographic_test_map.wcs.wcs.aux
    assert wcs_aux.hgln_obs is None
    assert wcs_aux.hglt_obs is None
    assert wcs_aux.dsun_obs is None


def test_rsun_meters_no_warning_for_hgs(heliographic_test_map):
    # Make sure that Stonyhurst heliographic maps do not emit a warning about assuming an
    # Earth-based observer when returning the physical radius of the Sun, because such an
    # assumption is not necessary

    # Convert the heliographic test map to Stonyhurst heliographic coordinates
    heliographic_test_map.meta.pop('HGLN_OBS')
    heliographic_test_map.meta.pop('HGLT_OBS')
    heliographic_test_map.meta.pop('DSUN_OBS')
    heliographic_test_map.meta['CTYPE1'] = 'HGLN-CAR'
    heliographic_test_map.meta['CTYPE2'] = 'HGLT-CAR'

    # Add a custom physical radius for the Sun
    heliographic_test_map.meta['rsun_ref'] = 1.1 * sunpy.sun.constants.radius.to_value(u.m)

    assert_quantity_allclose(heliographic_test_map.rsun_meters,
                             heliographic_test_map.meta['rsun_ref'] << u.m)


@figure_test
def test_rotation_rect_pixelated_data(aia171_test_map):
    aia_map = sunpy.map.Map(aia171_test_map)
    rect_map = aia_map.superpixel([2, 1] * u.pix, func=np.mean)
    rect_rot_map = rect_map.rotate(30 * u.deg)
    rect_rot_map.peek()


@pytest.mark.remote_data
@figure_test
def test_draw_contours_with_transform(sample_171, sample_hmi):
    aia_map = sunpy.map.Map(sample_171)
    hmi_map = sunpy.map.Map(sample_hmi)
    fig = plt.figure(figsize=(16, 4))

    # Panel 1: Implicit transform
    ax1 = fig.add_subplot(1, 3, 1, projection=aia_map)
    aia_map.plot(axes=ax1, clip_interval=(1, 99.99)*u.percent)
    hmi_map.draw_contours([-10, 10]*u.percent)
    ax1.set_title('Default, correct behavior')

    # Panel 2: Explicit transform
    ax2 = fig.add_subplot(1, 3, 2, projection=aia_map)
    aia_map.plot(axes=ax2, clip_interval=(1, 99.99)*u.percent)
    hmi_map.draw_contours([-10, 10]*u.percent, transform=ax2.get_transform(hmi_map.wcs))
    ax2.set_title('Explicitly specifying the correct transform')

    # Panel 3: Explicit transform with wacky rotation
    ax3 = fig.add_subplot(1, 3, 3, projection=aia_map)
    rotate_transform = Affine2D().rotate_deg_around(512, 512, 90)
    composite_transform = rotate_transform + ax3.get_transform(hmi_map.wcs)
    aia_map.plot(axes=ax3, clip_interval=(1, 99.99)*u.percent)
    hmi_map.draw_contours([-10, 10]*u.percent, transform=composite_transform)
    ax3.set_title('Contours rotated by 90 deg CCW')

    return fig


def test_plot_composite_map_updated_args(simple_map):
    simple_map.plot_settings['cmap'] = 'viridis'
    simple_map.plot_settings['norm'] = 'linear'
    simple_map.plot_settings['origin'] = 'upper'
    simple_map.plot_settings['alpha'] = 0.7
    simple_map.plot_settings['zorder'] = 8
    contour_args = {'norm': 'log',
                    'cmap':'plasma'}
    updated_args = simple_map._update_contour_args(contour_args)
    # Since 'norm' and  'cmap' are explicitly provided in contour_args of draw_contours,
    # their contour_args values will be used instead of plot_settings value
    assert updated_args ==  {
        'alpha': 0.7,
        'cmap': 'plasma',
        'norm': 'log',
        'origin': 'upper',
        'zorder': 8
    }


@figure_test
def test_draw_simple_map(simple_map):
    fig = plt.figure(figsize=(6, 6))
    ax = fig.add_subplot(1, 1, 1, projection=simple_map)
    simple_map.plot(axes=ax)
    return fig


@figure_test
def test_draw_carrington_map(carrington_map):
    fig = plt.figure(figsize=(6, 6))
    ax = fig.add_subplot(1, 1, 1, projection=carrington_map)
    carrington_map.plot(axes=ax)
    return fig


@pytest.mark.parametrize('method', _rotation_registry.keys())
@figure_test
def test_derotating_nonpurerotation_pcij(aia171_test_map, method):
    # The following map has a a PCij matrix that is not a pure rotation
    weird_map = aia171_test_map.rotate(30*u.deg).superpixel([2, 1]*u.pix)

    # De-rotating the map by its PCij matrix should result in a normal-looking map
    derotated_map = weird_map.rotate(method=method)

    fig = Figure(figsize=(8, 4))

    ax1 = fig.add_subplot(121, projection=weird_map)
    weird_map.plot(axes=ax1, title='Map with a non-pure-rotation PCij matrix')

    ax2 = fig.add_subplot(122, projection=derotated_map)
    derotated_map.plot(axes=ax2, title=f'De-rotated map via {method}')
    ax2.set_aspect(derotated_map.scale[1] / derotated_map.scale[0])

    return fig


# This function is used in the arithmetic tests below
def check_arithmetic_value_and_units(map_new, data_expected):
    assert u.allclose(map_new.quantity, data_expected)
    assert map_new.unit.is_equivalent(data_expected.unit)


@pytest.mark.parametrize('value', [
    10 * u.DN,
    u.Quantity([10], u.DN),
    u.Quantity(np.random.rand(128), u.DN),
    u.Quantity(np.random.rand(128, 128), u.DN),
])
def test_map_arithmetic_addition_subtraction(aia171_test_map, value):
    new_map = aia171_test_map + value
    check_arithmetic_value_and_units(new_map, aia171_test_map.quantity + value)
    new_map = value + aia171_test_map
    check_arithmetic_value_and_units(new_map, value + aia171_test_map.quantity)
    new_map = aia171_test_map - value
    check_arithmetic_value_and_units(new_map, aia171_test_map.quantity - value)
    new_map = value - aia171_test_map
    check_arithmetic_value_and_units(new_map, value - aia171_test_map.quantity)


@pytest.mark.parametrize('value', [
    10 * u.s,
    u.Quantity([10], u.s),
    u.Quantity(np.random.rand(128), u.s),
    u.Quantity(np.random.rand(128, 128), u.s),
    10.0,
    np.random.rand(128),
    np.random.rand(128, 128),
])
def test_map_arithmetic_multiplication_division(aia171_test_map, value):
    new_map = aia171_test_map * value
    check_arithmetic_value_and_units(new_map, aia171_test_map.quantity * value)
    new_map = value * aia171_test_map
    check_arithmetic_value_and_units(new_map, value * aia171_test_map.quantity)
    new_map = aia171_test_map / value
    check_arithmetic_value_and_units(new_map, aia171_test_map.quantity / value)
    with pytest.warns(RuntimeWarning, match='divide by zero encountered in'):  # NOQA: PT031
        new_map = value / aia171_test_map
        check_arithmetic_value_and_units(new_map, value / aia171_test_map.quantity)


def test_map_arithmetic_pow(aia171_test_map):
    new_map = aia171_test_map ** 2
    check_arithmetic_value_and_units(new_map, aia171_test_map.quantity ** 2)


def test_map_arithmetic_neg(aia171_test_map):
    new_map = -aia171_test_map
    check_arithmetic_value_and_units(new_map, -aia171_test_map.quantity)


@pytest.mark.parametrize("value", ['map', 'foobar', None, ['foo', 'bar']])
def test_map_arithmetic_operations_raise_exceptions(aia171_test_map, value):
    value = aia171_test_map if value == 'map' else value
    with pytest.raises(TypeError):
        _ = aia171_test_map + value
    with pytest.raises(TypeError):
        _ = aia171_test_map * value
    with pytest.raises(TypeError):
        _ = value / aia171_test_map

@pytest.mark.parametrize(('units_string','expected_unit'),[
    ('Gauss', u.G),
    ('G', u.G),
    ('DN', u.DN),
    ('DN/s', u.DN/u.s),
    ('DN/pix', u.DN/u.pixel),
    ('DN / pix', u.DN/u.pixel),
    ('DN sr / s', u.DN*u.sr/u.s),
    ('DN/(pix s)', u.DN/u.pixel/u.s),
    ('counts / pixel', u.ct/u.pix),
])
def test_parse_fits_units(units_string, expected_unit):
    out_unit = GenericMap._parse_fits_unit(units_string)
    assert out_unit == expected_unit


@pytest.mark.parametrize('units_string', ['DN / electron', 'electron', 'Mx'])
def test_parse_nonfits_units(units_string):
    with pytest.warns(SunpyMetadataWarning, match='Could not parse unit string'):
        assert GenericMap._parse_fits_unit(units_string) is None


def test_only_cd():
    data = np.ones([6, 6], dtype=np.float64)
    header = {
        'CRVAL1': 0,
        'CRVAL2': 0,
        'CRPIX1': 5,
        'CRPIX2': 5,
        'CD1_1': 3,
        'CD1_2': -4,
        'CD2_1': 5,
        'CD2_2': 12,
        'NAXIS1': 6,
        'NAXIS2': 6,
        'CUNIT1': 'arcsec',
        'CUNIT2': 'arcsec',
        'CTYPE1': 'HPLN-TAN',
        'CTYPE2': 'HPLT-TAN',
    }
    cd_map = sunpy.map.Map((data, header))
    np.testing.assert_allclose(u.Quantity(cd_map.scale).value, np.array([5, 13]))
    np.testing.assert_allclose(cd_map.rotation_matrix, np.array([[3/5, -4/5], [5/13, 12/13]]))


def test_submap_nan_error(aia171_test_map):
    # See https://github.com/sunpy/sunpy/pull/7543#issuecomment-2167019208 for more context
    coord_native = SkyCoord(0*u.arcsec, 0*u.arcsec, frame=aia171_test_map.coordinate_frame)
    aia171_test_map.submap(coord_native, width=1000*u.arcsec, height=1000*u.arcsec)
    coord_other = SkyCoord(0*u.arcsec, 0*u.arcsec, frame='helioprojective', observer='earth', obstime=aia171_test_map.date)
    with pytest.raises(ValueError, match="The provided input coordinates to"):
        aia171_test_map.submap(coord_other, width=1000*u.arcsec, height=1000*u.arcsec)
