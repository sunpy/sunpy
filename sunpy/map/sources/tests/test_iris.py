import pytest

import astropy.units as u

from sunpy.data.test import get_dummy_map_from_header, get_test_filepath
from sunpy.map.sources.iris import SJIMap
from sunpy.util.exceptions import SunpyMetadataWarning

__author__ = 'Pritish C. (VaticanCameos)'


@pytest.fixture
def iris_map():
    header = get_test_filepath("iris_l2_20130801_074720_4040000014_SJI_1400_t000.header")
    return get_dummy_map_from_header(header)


def test_fitstoIRIS(iris_map):
    assert isinstance(iris_map, SJIMap)


def test_is_datasource_for(iris_map):
    assert iris_map.is_datasource_for(iris_map.data, iris_map.meta)


def test_observatory(iris_map):
    assert iris_map.observatory == "IRIS"


def test_wavelength(iris_map):
    assert iris_map.wavelength == u.Quantity(1400, 'Angstrom')


def test_level_number(iris_map):
    assert iris_map.processing_level == 2.0


def test_units(iris_map):
    assert iris_map.unit == u.ct


def test_wcs(iris_map):
    # Smoke test that WCS is valid and can transform from pixels to world coordinates
    with pytest.warns(SunpyMetadataWarning, match='Missing metadata for observer'):
        iris_map.pixel_to_world(0*u.pix, 0*u.pix)
