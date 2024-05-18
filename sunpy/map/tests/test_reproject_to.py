"""
Test the `GenericMap.reproject_to()` method
"""
import warnings

import numpy as np
import pytest
from matplotlib.testing.decorators import check_figures_equal

import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.wcs import WCS

import sunpy.coordinates
import sunpy.map
from sunpy.coordinates._transformations import propagate_with_solar_surface
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
    with pytest.warns(SunpyUserWarning, match="rsun mismatch detected: "):
        # Modifying the `hpc_header` rsun value to create a mismatch
        hpc_header["rsun_ref"] += 1

        # Reproject with the mismatched rsun
        aia171_test_map.reproject_to(hpc_header)


def test_reproject_to_warn_using_contexts(aia171_test_map, hpc_header):
    with propagate_with_solar_surface():
        with sunpy.coordinates.SphericalScreen(aia171_test_map.observer_coordinate):
            # Check if a warning is raised if both context managers are used at the same time.
            with pytest.warns(UserWarning, match="Using propagate_with_solar_surface and SphericalScreen together result in loss of off-disk data."):
                aia171_test_map.reproject_to(hpc_header)
