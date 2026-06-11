"""
Test suite for METISMap class.

This module provides unit tests for the METISMap source, ensuring correct
metadata handling, coordinate frame definition, and instrument-specific
functionalities like occultor masking.
"""

import numpy as np
import pytest

import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.visualization import ImageNormalize

from sunpy.map import Map
from sunpy.map.sources.solo import METISMap
from sunpy.util.exceptions import SunpyUserWarning

METIS_HEADER_VARIANTS = [
    pytest.param({}, id="baseline-VL"),
    pytest.param({"FILTER": "UV", "BTYPE": "UV Lyman-alpha intensity"}, id="UV-filter"),
    pytest.param({"HGLN_OBS": 45.0, "HGLT_OBS": -7.25}, id="off-axis-observer"),
    pytest.param({"DSUN_OBS": 4.5e10}, id="perihelion"),
    pytest.param({"CRVAL1": 500.0, "CRVAL2": 300.0}, id="off-centre-pointing"),
]

_BASE_METIS_HEADER = {
    "INSTRUME": "Metis",
    "OBSRVTRY": "Solar Orbiter",
    "FILTER": "VL",
    "BTYPE": "VL total brightness",
    "LEVEL": "L2",
    "DATE-AVG": "2023-01-01T12:00:00",
    "CDELT1": 10.0,
    "CDELT2": 10.0,
    "CUNIT1": "arcsec",
    "CUNIT2": "arcsec",
    "CTYPE1": "HPLN-TAN",
    "CTYPE2": "HPLT-TAN",
    "NAXIS1": 1024,
    "NAXIS2": 1024,
    "CRPIX1": 512.0,
    "CRPIX2": 512.0,
    "CRVAL1": 0.0,
    "CRVAL2": 0.0,
    "RSUN_ARC": 960.0,
    "DSUN_OBS": 1.5e11,
    "HGLN_OBS": 0.0,
    "HGLT_OBS": 0.0,
    "INN_FOV": 1.6,
    "OUT_FOV": 2.9,
    "IO_XCEN": 512.0,
    "IO_YCEN": 512.0,
    "FS_XCEN": 512.0,
    "FS_YCEN": 512.0,
    "SUN_XCEN": 512.0,
    "SUN_YCEN": 512.0,
}


@pytest.fixture(scope="module", params=METIS_HEADER_VARIANTS)
def metis_map(request,metis_test_data):
    header = {**_BASE_METIS_HEADER, **request.param}
    return Map(metis_test_data, header)


@pytest.fixture
def minimal_metis_header():
    return _BASE_METIS_HEADER.copy()


@pytest.fixture
def metis_map_non_square(metis_test_data, minimal_metis_header):
    header = {**minimal_metis_header, "CDELT1": 10.0, "CDELT2": 20.0}
    return Map(metis_test_data, header)


@pytest.fixture(scope="module")
def metis_test_data():
    return np.zeros((1024, 1024), dtype=np.float32)

# basic tests
def test_basic_initialization(metis_map):
    assert isinstance(metis_map, METISMap)
    assert metis_map.instrument == "Metis"


def test_observatory(metis_map):
    assert metis_map.observatory == "Solar Orbiter"


def test_coordinate_frame(metis_map):
    assert metis_map.coordinate_frame.name == "helioprojective"


def test_is_datasource_for(metis_map):
    assert metis_map.is_datasource_for(metis_map.data, metis_map.meta)

# wcs
def test_wcs(metis_map):
    metis_map.pixel_to_world(0 * u.pix, 0 * u.pix)


def test_wcs_center_pixel(metis_test_data, minimal_metis_header):
    """Test that pixel coordinates correctly transform to world coordinates (Arcsecs)."""
    metis_map = Map(metis_test_data,minimal_metis_header)
    # Center of the image (0-indexed 511.5) should be close to CRVAL (0,0)
    center_coord = metis_map.pixel_to_world(511 * u.pix, 511 * u.pix)
    assert isinstance(center_coord, SkyCoord)
    assert u.allclose(center_coord.Tx, 0 * u.arcsec, atol=1e-2 * u.arcsec)


# plotting
def test_cmap_by_filter(metis_map):
    """Verify that the default colormap is assigned based on the FILTER keyword."""
    expected = {
        "VL": "solometisvl-tb",
        "UV": "solometisuv",
    }[metis_map.meta["FILTER"]]
    assert metis_map.plot_settings["cmap"] == expected


def test_norm(metis_map):
    assert isinstance(metis_map.plot_settings["norm"], ImageNormalize)

# measurements
@pytest.mark.parametrize(("btype", "expected"), [
    ("VL total brightness",             "VL-TB"),
    ("VL polarized brightness",         "VL-PB"),
    ("VL fixed-polarization intensity", "VL-FP"),
    ("VL polarization angle",           "VL-PA"),
    ("Stokes I",                        "VL-SI"),
    ("Stokes Q",                        "VL-SQ"),
    ("Stokes U",                        "VL-SU"),
    ("UV Lyman-alpha intensity",        "UV"),
])
def test_measurement(metis_test_data, minimal_metis_header, btype, expected):
    filter_key = "UV" if expected == "UV" else "VL"
    header = {**minimal_metis_header, "FILTER": filter_key, "BTYPE": btype}
    m = Map(metis_test_data, header)
    assert m.measurement == expected


# rsun things
def test_rsun_obs_uses_super(metis_test_data, minimal_metis_header):
    """Verify rsun_obs uses GenericMap when RSUN_OBS exists."""
    header = {**minimal_metis_header, "RSUN_OBS": 950.0}
    metis_map = Map(metis_test_data, header)
    assert metis_map.rsun_obs == 950.0 * u.arcsec


def test_fall_back_to_rsun_arc(metis_test_data, minimal_metis_header):
    """Verify rsun_obs falls back to RSUN_ARC when no other solar radius keywords exist."""
    metis_map = Map(metis_test_data, minimal_metis_header)
    assert metis_map.rsun_obs == 960.0 * u.arcsec

# testing the mask
def test_mask_is_set(metis_map):
    assert metis_map.mask is not None


def test_mask_dtype(metis_map):
    """Verify mask was created."""
    assert metis_map.mask.dtype == bool


def test_mask_occ_does_not_mask_entire_image(metis_map):
    """Verify that mask_occ does not occult the entire image"""
    assert not metis_map.mask.all()


def test_mask_occ_center_occulted(metis_map):
    """Verify that the center is masked"""
    assert metis_map.mask[512,512]


def test_mask_non_square_pixels_returns_none(metis_map_non_square):
    result = metis_map_non_square.mask
    assert result is None


def test_mask_missing_keys_warns_and_returns_none(metis_test_data, minimal_metis_header):
    """Verify that mask falls back to super().mask when occulter keys are missing."""
    header = minimal_metis_header.copy()
    del header["INN_FOV"]  # Remove one required key to trigger missing_keys
    metis_map = Map(metis_test_data, header)
    with pytest.warns(SunpyUserWarning, match="Missing.*keys required to calculate occulter mask"):
        result = metis_map.mask
    # GenericMap.mask returns None by default when no mask is set
    assert result is None


def test_mask_dr1_workaround(metis_test_data, minimal_metis_header):
    # data release 1: fs_xcen == crpix1, so sun_xcen/ycen should be used
    header = {**minimal_metis_header, "FS_XCEN": 512.0, "FS_YCEN": 512.0,
              "CRPIX1": 512.0, "CRPIX2": 512.0, "SUN_XCEN": 512.0, "SUN_YCEN": 512.0}
    m = Map(metis_test_data, header)
    assert m.mask is not None
    assert m.mask.dtype == bool
    assert not m.mask.all()


def test_mask_uses_field_stop_center_dr2(metis_test_data, minimal_metis_header):
    """Verify the normal (data release 2) mask path where fs_xcen != crpix1 or fs_ycen != crpix2."""
    header = {**minimal_metis_header, "FS_XCEN": 514.0, "FS_YCEN": 514.0}
    metis_map = Map(metis_test_data, header)
    assert metis_map.mask is not None
    assert metis_map.mask.dtype == bool
    assert not metis_map.mask.all()
