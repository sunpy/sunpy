"""
Test the `GenericMap.reproject_to()` method
"""
import warnings

import numpy as np
import pytest

import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.wcs import WCS

import sunpy.map
from sunpy.tests.helpers import figure_test


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
@pytest.mark.filterwarnings('ignore:Missing metadata')  # TODO: fix bug for HGS maps
def test_reproject_to_hgs(aia171_test_map, hgs_header):
    aia171_test_map.reproject_to(hgs_header).plot()


@figure_test
@pytest.mark.filterwarnings('ignore:Missing metadata')  # TODO: fix bug for HGS maps
def test_reproject_to_hgs_wcs(aia171_test_map, hgs_header):
    aia171_test_map.reproject_to(WCS(hgs_header)).plot()


@figure_test
def test_reproject_to_hpc(aia171_test_map, hpc_header):
    aia171_test_map.reproject_to(hpc_header).plot()


@figure_test
def test_reproject_to_hpc_interpolation(aia171_test_map, hpc_header):
    aia171_test_map.reproject_to(hpc_header, algorithm='interpolation').plot()


@figure_test
def test_reproject_to_hpc_exact(aia171_test_map, hpc_header):
    aia171_test_map.reproject_to(hpc_header, algorithm='exact').plot()


@figure_test
def test_reproject_to_hpc_adaptive(aia171_test_map, hpc_header):
    aia171_test_map.reproject_to(hpc_header, algorithm='adaptive').plot()


def test_return_footprint(aia171_test_map, hpc_header):
    pytest.importorskip("reproject")

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
    pytest.importorskip("reproject")
    with pytest.raises(ValueError):
        aia171_test_map.reproject_to(hpc_header, algorithm='something')
