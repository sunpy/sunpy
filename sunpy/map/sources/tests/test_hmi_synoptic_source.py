import pytest

import astropy.units as u

from sunpy.data.test import get_dummy_map_from_header, get_test_filepath
from sunpy.map.sources.sdo import HMISynopticMap
from sunpy.util.exceptions import SunpyMetadataWarning
from .helpers import _test_private_date_setters


@pytest.fixture
def hmi_synoptic():
    return get_dummy_map_from_header(get_test_filepath('hmi_synoptic.header'))


def test_fitstoHMISynoptic(hmi_synoptic):
    """Tests the creation of HMISynopticMap using FITS."""
    assert isinstance(hmi_synoptic, HMISynopticMap)


def test_is_datasource_for(hmi_synoptic):
    """Test the is_datasource_for method of HMISynopticMap.
    Note that header data to be provided as an argument
    can be a MetaDict object, which in this case is
    hmi.meta."""
    assert hmi_synoptic.is_datasource_for(hmi_synoptic.data, hmi_synoptic.meta)


def test_observatory(hmi_synoptic):
    """Tests the observatory property of the HMISynopticMap object."""
    assert hmi_synoptic.observatory == "SDO"


def test_measurement(hmi_synoptic):
    """Tests the measurement property of the HMISynopticMap object."""
    assert hmi_synoptic.measurement == "carrington"


def test_reference_date(hmi_synoptic):
    assert hmi_synoptic.reference_date.isot == "2018-11-09T12:30:52.000"


def test_date(hmi_synoptic):
    assert hmi_synoptic.date.isot == "2018-11-09T12:30:52.000"


def test_private_date_setters(hmi_synoptic):
    _test_private_date_setters(hmi_synoptic)


def test_unit(hmi_synoptic):
    assert hmi_synoptic.unit == u.G
    assert hmi_synoptic.unit == u.Unit("Mx/cm^2")
    assert hmi_synoptic.unit.to_string() == 'Mx / cm2'
    hmi_synoptic.meta['bunit'] = 'm'
    assert hmi_synoptic.unit == u.m


def test_wcs(hmi_synoptic):
    # Smoke test that WCS is valid and can transform from pixels to world coordinates
    with pytest.warns(SunpyMetadataWarning, match='Missing metadata for observer'):
        hmi_synoptic.pixel_to_world(0*u.pix, 0*u.pix)
