
import pytest

import astropy.units as u

from sunpy.data.test import get_dummy_map_from_header, get_test_filepath
from sunpy.map.mapbase import SpatialPair
from sunpy.map.sources.hinode import SOTMap
from sunpy.util.exceptions import SunpyMetadataWarning
from .helpers import _test_private_date_setters


@pytest.fixture
def sot():
    return get_dummy_map_from_header(get_test_filepath("HinodeSOT.header"))


def test_fitstoSOT(sot):
    """Tests the creation of SOTMap using FITS."""
    assert isinstance(sot, SOTMap)


def test_sot_coordinate_system(sot):
    assert sot.coordinate_system ==  SpatialPair(axis1='HPLN-TAN', axis2='HPLT-TAN')


def test_is_datasource_for(sot):
    """Test the is_datasource_for method of SOTMap.
    Note that header data to be provided as an argument
    can be a MetaDict object."""
    assert sot.is_datasource_for(sot.data, sot.meta)


def test_observatory(sot):
    """Tests the observatory property of the SOTMap object."""
    assert sot.observatory == "Hinode"


def test_detector(sot):
    """Tests the detector property of the SOTMap object"""
    assert sot.detector == "SOT"


def test_reference_date(sot):
    assert sot.reference_date.isot == "2015-10-13T23:13:44.601"


def test_date(sot):
    assert sot.date.isot == "2015-10-13T23:13:44.601"


def test_private_date_setters(sot):
    _test_private_date_setters(sot)


def test_measurement(sot):
    """Tests the measurement property of the SOTMap object."""
    assert sot.measurement is None


def test_instruments(sot):
    """Tests the Instruments object of SOTMap."""
    assert (sot.Instruments == ['SOT/WB',
                                'SOT/NB', 'SOT/SP', 'SOT/CT'])


def test_waves(sot):
    """Tests the Waves object of SOTMap."""
    assert (sot.Waves == ['6302A', 'BFI no move',
                          'CN bandhead 3883', 'Ca II H line',
                          'G band 4305', 'NFI no move', 'TF Fe I 6302',
                          'TF Mg I 5172', 'TF Na I 5896',
                          'blue cont 4504', 'green cont 5550',
                          'red cont 6684'])


def test_obstype(sot):
    """Tests the Observation_Type object of SOTMap."""
    assert (sot.Observation_Type == ['FG (simple)',
                                     'FG focus scan', 'FG shuttered I and V',
                                     'FG shutterless I and V', 'FG shutterless I and V with 0.2s intervals',
                                     'FG shutterless Stokes', 'SP IQUV 4D array'])


def test_wcs(sot):
    # Smoke test that WCS is valid and can transform from pixels to world coordinates
    with pytest.warns(SunpyMetadataWarning, match='assuming Earth-based observer'):
        sot.pixel_to_world(0*u.pix, 0*u.pix)
