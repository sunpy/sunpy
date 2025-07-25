"""
Test the `GenericMap.reproject_to()` method
"""
import warnings

import numpy as np
import pytest
from matplotlib.figure import Figure
from matplotlib.testing.decorators import check_figures_equal

import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.tests.helper import assert_quantity_allclose
from astropy.wcs import WCS

import sunpy.coordinates
import sunpy.map
from sunpy.tests.helpers import figure_test
from sunpy.util.exceptions import SunpyUserWarning


@pytest.fixture
def hgs_header(aia171_test_map):
    return sunpy.map.make_fitswcs_header(
        (180, 360),
        SkyCoord(0*u.deg, 0*u.deg,
                 frame='heliographic_stonyhurst',
                 obstime=aia171_test_map.date,
                 rsun=aia171_test_map.coordinate_frame.rsun),
        scale=(1, 1)*u.deg/u.pix,
        projection_code='CAR'
    )


@pytest.fixture
def hpc_header(aia171_test_map):
    new_observer = SkyCoord(45*u.deg, 0*u.deg, 1*u.AU,
                            frame='heliographic_stonyhurst',
                            obstime=aia171_test_map.date)
    return sunpy.map.make_fitswcs_header(
        aia171_test_map.data.shape,
        SkyCoord(0*u.arcsec, 0*u.arcsec,
                 frame='helioprojective',
                 obstime=aia171_test_map.date,
                 observer=new_observer,
                 rsun=aia171_test_map.coordinate_frame.rsun),
        scale=u.Quantity(aia171_test_map.scale),
        projection_code='TAN'
    )


@figure_test
def test_reproject_to_hgs(aia171_test_map, hgs_header):
    aia171_test_map.reproject_to(hgs_header).plot()


@check_figures_equal(extensions=["png"])
def test_reproject_to_hgs_wcs(fig_test, fig_ref, aia171_test_map, hgs_header):
    with warnings.catch_warnings():
        # NumPy <1.19 emits a RuntimeWarning because of comparison against NaNs
        warnings.filterwarnings("ignore", message='invalid value encountered',
                                category=RuntimeWarning)

        # Tests whether reprojecting to a WCS instance gives the same answer as to a header
        header_map = aia171_test_map.reproject_to(hgs_header)
        wcs_map = aia171_test_map.reproject_to(WCS(hgs_header))

        ax_ref = fig_ref.add_subplot(projection=header_map)
        header_map.plot(axes=ax_ref)

        ax_test = fig_test.add_subplot(projection=wcs_map)
        wcs_map.plot(axes=ax_test)


@check_figures_equal(extensions=["png"])
def test_reproject_to_hpc_default(fig_test, fig_ref, aia171_test_map, hpc_header):
    with warnings.catch_warnings():
        # NumPy <1.19 emits a RuntimeWarning because of comparison against NaNs
        warnings.filterwarnings("ignore", message='invalid value encountered',
                                category=RuntimeWarning)

        # Tests whether the default reprojection is "interpolation"
        default_map = aia171_test_map.reproject_to(hpc_header)
        interpolation_map = aia171_test_map.reproject_to(hpc_header, algorithm='interpolation')

        ax_ref = fig_ref.add_subplot(projection=interpolation_map)
        interpolation_map.plot(axes=ax_ref)

        ax_test = fig_test.add_subplot(projection=default_map)
        default_map.plot(axes=ax_test)


@figure_test
def test_reproject_to_hpc_interpolation(aia171_test_map, hpc_header):
    aia171_test_map.reproject_to(hpc_header, algorithm='interpolation').plot()


@figure_test
def test_reproject_to_hpc_exact(aia171_test_map, hpc_header):
    aia171_test_map.reproject_to(hpc_header, algorithm='exact').plot()


@figure_test
def test_reproject_to_hpc_adaptive(aia171_test_map, hpc_header):
    aia171_test_map.reproject_to(hpc_header, algorithm='adaptive', kernel='Hann', boundary_mode='strict').plot()


def test_return_footprint(aia171_test_map, hpc_header):
    with warnings.catch_warnings():
        # NumPy <1.19 emits a RuntimeWarning because of comparison against NaNs
        warnings.filterwarnings("ignore", message='invalid value encountered',
                                category=RuntimeWarning)

        return_without_footprint = aia171_test_map.reproject_to(hpc_header)
        assert isinstance(return_without_footprint, sunpy.map.GenericMap)

        return_with_footprint = aia171_test_map.reproject_to(hpc_header, return_footprint=True)
        assert len(return_with_footprint) == 2
        assert isinstance(return_with_footprint[0], sunpy.map.GenericMap)
        assert isinstance(return_with_footprint[1], np.ndarray)


def test_invalid_inputs(aia171_test_map, hpc_header):
    with pytest.raises(ValueError, match="The specified algorithm must be one of"):
        aia171_test_map.reproject_to(hpc_header, algorithm='something')


def test_rsun_mismatch_warning(aia171_test_map, hpc_header):
    # Modifying the `hpc_header` rsun value to create a mismatch
    hpc_header["rsun_ref"] += 1
    with pytest.warns(SunpyUserWarning, match="rsun mismatch detected: "):
        # Reproject with the mismatched rsun
        aia171_test_map.reproject_to(hpc_header)


@figure_test
@pytest.mark.parametrize('auto_extent', [None, 'corners', 'edges', 'all'])
def test_reproject_to_auto_extent(aia171_test_map, auto_extent):
    aia = aia171_test_map.submap([32, 32] * u.pix, top_right=[80, 100]*u.pix)

    header = {
        'naxis1': 60,
        'naxis2': 200,
        'crpix1': 30.5,
        'crpix2': 100.5,
        'crval1': -750,
        'crval2': 0,
        'cdelt1': 10,
        'cdelt2': 10,
        'cunit1': 'arcsec',
        'cunit2': 'arcsec',
        'ctype1': 'HPLN-TAN',
        'ctype2': 'HPLT-TAN',
        'hgln_obs': 80,
        'hglt_obs': 10,
        'rsun_ref': aia.rsun_meters.value,
        'dsun_obs': aia.dsun.value,
        'date-obs': aia.date.isot,
        'mjd-obs': aia.date.mjd,
    }

    reprojected_aia = aia.reproject_to(header, auto_extent=auto_extent)
    fig = Figure()
    ax = fig.add_subplot(projection=reprojected_aia)
    reprojected_aia.plot(axes=ax, title=f"auto_extent={auto_extent}")
    aia.draw_extent(axes=ax, color='red')
    return fig


@pytest.mark.parametrize(('auto_extent', 'crpix', 'pixel_shape'), [(None, [5.5, 5.5], [10, 10]),
                                                                   ('corners', [0.5, 836.5], [241, 1672]),
                                                                   ('edges', [0.5, 836.5], [482, 1672]),
                                                                   ('all', [0.5, 836.5], [963, 1672])])
def test_reproject_to_auto_extent_wcs(aia171_test_map, auto_extent, crpix, pixel_shape):
    hgs_header = {
        'naxis1': 15,
        'naxis2': 12,
        'crpix1': 0.5,
        'crpix2': 6.5,
        'crval1': 0,
        'crval2': 0,
        'cdelt1': 10,
        'cdelt2': 10,
        'cunit1': 'deg',
        'cunit2': 'deg',
        'ctype1': 'HGLN-CAR',
        'ctype2': 'HGLT-CAR',
        'rsun_ref': 700000000,
        'dsun_obs': 150000000000,
        'date-obs': '2001-01-01',
    }
    hgs_map = sunpy.map.Map((np.arange(180).reshape((12, 15)), hgs_header))

    hpc_header = {
        'naxis1': 10,
        'naxis2': 10,
        'crpix1': 5.5,
        'crpix2': 5.5,
        'crval1': 0,
        'crval2': 0,
        'cdelt1': 1,
        'cdelt2': 1,
        'cunit1': 'arcsec',
        'cunit2': 'arcsec',
        'ctype1': 'HPLN-TAN',
        'ctype2': 'HPLT-TAN',
        'hgln_obs': 0,
        'hglt_obs': 0,
        'rsun_ref': 700000000,
        'dsun_obs': 150000000000,
        'date-obs': '2001-01-01',
        'mjd-obs': 51910,
    }

    reprojected_wcs = hgs_map.reproject_to(hpc_header, auto_extent=auto_extent).wcs

    # These parts of the WCS should not change
    assert_quantity_allclose(reprojected_wcs.wcs.crval, hgs_map.wcs.wcs.crval)
    assert_quantity_allclose(reprojected_wcs.wcs.cdelt * u.deg, 1 * u.arcsec)

    # These parts of the WCS should have changed
    assert_quantity_allclose(reprojected_wcs.wcs.crpix, crpix)
    assert_quantity_allclose(reprojected_wcs.pixel_shape, pixel_shape)


@figure_test
@pytest.mark.parametrize("screen", [sunpy.coordinates.SphericalScreen, sunpy.coordinates.PlanarScreen])
def test_reproject_to_screen_plus_diffrot(aia171_test_map, screen):
    new_time = aia171_test_map.date + 5*u.day
    header = {
        'naxis1': 150,
        'naxis2': 150,
        'crpix1': 75.5,
        'crpix2': 75.5,
        'crval1': 0,
        'crval2': 0,
        'cdelt1': 20,
        'cdelt2': 20,
        'cunit1': 'arcsec',
        'cunit2': 'arcsec',
        'ctype1': 'HPLN-TAN',
        'ctype2': 'HPLT-TAN',
        'hgln_obs': 0,
        'hglt_obs': 30,
        'rsun_ref': aia171_test_map.rsun_meters.value,
        'dsun_obs': aia171_test_map.dsun.value,
        'date-obs': new_time.isot,
        'mjd-obs': new_time.mjd,
    }

    with sunpy.coordinates.transform_with_sun_center(), screen(aia171_test_map.observer_coordinate, only_off_disk=True):
        without_diffrot = aia171_test_map.reproject_to(header)
        with sunpy.coordinates.propagate_with_solar_surface():
            with_diffrot = aia171_test_map.reproject_to(header)

    fig = Figure(figsize=(11, 4))
    ax1 = fig.add_subplot(121, projection=without_diffrot)
    without_diffrot.plot(axes=ax1, title="Without diffrot")
    ax2 = fig.add_subplot(122, projection=with_diffrot)
    with_diffrot.plot(axes=ax2, title="With diffrot")
    return fig
