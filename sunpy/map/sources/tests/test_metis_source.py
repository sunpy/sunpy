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

from sunpy.map import Map
from sunpy.map.sources.solo import METISMap

# ============================================================================
# FIXTURES - Test data and metadata configuration
# ============================================================================

METIS_HEADER_VARIANTS = [
    pytest.param({}, id="baseline-VL"),
    pytest.param({"FILTER": "UV", "BTYPE": "UV Lyman-alpha intensity"}, id="UV-filter"),
    pytest.param({"HGLN_OBS": 45.0, "HGLT_OBS": -7.25}, id="off-axis-observer"),
    pytest.param({"DSUN_OBS": 4.5e10}, id="perihelion"),
    pytest.param({"CRVAL1": 500.0, "CRVAL2": 300.0}, id="off-centre-pointing"),
]

_BASE_METIS_HEADER = {
    "INSTRUME": "METIS",
    "OBSRVTRY": "SOLAR ORBITER",
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


@pytest.fixture(scope="module")
def metis_test_data():
    """
    Generates synthetic 1024x1024 image data with a radial gradient.
    """
    size = 1024
    y, x = np.ogrid[-size / 2 : size / 2, -size / 2 : size / 2]
    r = np.sqrt(x**2 + y**2)
    # Exponential decay to simulate coronal brightness
    data = 1000 * np.exp(-r / 200) + np.random.normal(10, 2, (size, size))
    return data.astype(np.float32)


def test_coordinate_frame( metis_map):
    """
    Verify that the map is correctly assigned a Helioprojective frame
    and the observer is identified as Solar Orbiter.
    """
    assert metis_map.coordinate_frame.name == "helioprojective"
    assert metis_map.observatory == "SOLAR ORBITER"


def test_colormaps_by_filter( metis_map):
    """
    Verify that the default colormap is assigned based on the FILTER keyword.
    """
    expected = {
        "VL": "solometisvl-tb",
        "UV": "solometisuv",
    }[metis_map.meta["FILTER"]]
    assert metis_map.plot_settings["cmap"] == expected


def test_wcs_conversion( metis_test_data, minimal_metis_header):
    """
    Test that pixel coordinates correctly transform to world coordinates (Arcsecs).
    """
    metis_map = Map(metis_test_data,minimal_metis_header)
    # Center of the image (0-indexed 511.5) should be close to CRVAL (0,0)
    center_coord = metis_map.pixel_to_world(511 * u.pix, 511 * u.pix)
    assert isinstance(center_coord, SkyCoord)
    assert u.allclose(center_coord.Tx, 0 * u.arcsec, atol=1e-2 * u.arcsec)


"""
Tests for proper object instantiation.
"""

def test_basic_initialization( metis_map):
    assert isinstance(metis_map, METISMap)
    assert metis_map.instrument == "METIS"


def test_rsun_fallback_logic( metis_test_data, minimal_metis_header):
    """
    Verify the hierarchical search for solar radius keywords
    (RSUN_ARC vs RSUN_OBS).
    """
    header = minimal_metis_header.copy()
    del header["RSUN_ARC"]
    header["RSUN_OBS"] = 950.0
    metis_map = Map(metis_test_data, header)
    # SunPy normalizes this to rsun_obs in metadata
    assert metis_map.meta["rsun_obs"] == 950.0


def test_prodtype_property(metis_test_data, minimal_metis_header):
    """Verify that prodtype is correctly derived from BTYPE."""
    metis_map = Map(metis_test_data,minimal_metis_header)
    assert metis_map.measurement == "VL-TB"


def test_rsun_obs_uses_super(metis_test_data, minimal_metis_header):
    """
    Verify rsun_obs uses GenericMap when RSUN_OBS exists.
    """
    header = minimal_metis_header.copy()
    header["RSUN_OBS"] = 950.0

    metis_map = Map(metis_test_data, header)

    assert metis_map.rsun_obs == 950.0 * u.arcsec

def test_fall_back_to_rsun_arc(metis_test_data, minimal_metis_header):
    """
    Verify rsun_obs falls back to RSUN_ARC when no other solar radius keywords exist.
    """
    header = minimal_metis_header.copy()

    header.pop("RSUN_OBS", None)
    header.pop("SOLAR_R", None)
    header.pop("RADIUS", None)

    metis_map = Map(metis_test_data, header)

    assert metis_map.rsun_obs == 960.0 * u.arcsec
