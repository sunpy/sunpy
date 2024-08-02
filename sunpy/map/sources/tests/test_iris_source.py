"""
Test cases for IRIS SJIMap
"""
import pytest

import astropy.units as u

from sunpy.data.test import get_dummy_map_from_header, get_test_filepath
from sunpy.map.sources.iris import SJIMap
from sunpy.util.exceptions import SunpyMetadataWarning
from .helpers import _test_private_date_setters

__author__ = 'Pritish C. (VaticanCameos)'


@pytest.fixture
def irismap():
    header = get_test_filepath("iris_l2_20130801_074720_4040000014_SJI_1400_t000.header")
    return get_dummy_map_from_header(header)


def test_fitstoIRIS(irismap):
    """Tests the creation of SJIMap using FITS."""
    assert (isinstance(irismap, SJIMap))


def test_private_date_setters(irismap):
    _test_private_date_setters(irismap)


def test_is_datasource_for(irismap):
    """Test the is_datasource_for method of SJIMap.
    Note that header data to be provided as an argument
    can be a MetaDict object."""
    assert irismap.is_datasource_for(irismap.data, irismap.meta)


def test_observatory(irismap):
    """Tests the observatory property of SJIMap."""
    assert irismap.observatory == "IRIS"


def test_wavelength(irismap):
    """Tests the wavelength and waveunit property of the SJIMap"""
    assert irismap.wavelength == u.Quantity(1400, 'Angstrom')


def test_level_number(irismap):
    """Tests the processing_level property of the SJIMap"""
    assert irismap.processing_level == 2.0


def test_units(irismap):
    """Tests the unit property of the SJIMap"""
    assert irismap.unit == u.Unit("DN")


def test_wcs(irismap):
    # Smoke test that WCS is valid and can transform from pixels to world coordinates
    with pytest.warns(SunpyMetadataWarning, match='Missing metadata for observer'):
        irismap.pixel_to_world(0*u.pix, 0*u.pix)
