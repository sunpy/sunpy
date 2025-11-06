"""
Test cases for SOHO MDIMap subclass.
"""
import pytest

import astropy.units as u
from astropy.coordinates import Angle

from sunpy.data.test import get_dummy_map_from_header, get_test_filepath
from sunpy.map.sources.soho import MDIMap, MDISynopticMap
from sunpy.util.exceptions import SunpyMetadataWarning
from .helpers import _test_private_date_setters

__author__ = 'Pritish C. (VaticanCameos)'


@pytest.fixture
def mdi():
    return get_dummy_map_from_header(get_test_filepath("mdi.fd_Ic.20101015_230100_TAI.data.header"))


@pytest.fixture
def mdi_synoptic():
    return get_dummy_map_from_header(get_test_filepath('mdi_synoptic.header'))


def test_fitstoMDI(mdi):
    """Tests the creation of MDIMap using FITS."""
    assert isinstance(mdi, MDIMap)


def test_is_datasource_for(mdi):
    """Test the is_datasource_for method of MDIMap.
    Note that header data to be provided as an argument
    can be a MetaDict object."""
    assert mdi.is_datasource_for(mdi.data, mdi.meta)


def test_observatory(mdi):
    """Tests the observatory property of the MDIMap object."""
    assert mdi.observatory == "SOHO"


def test_reference_date(mdi):
    assert mdi.reference_date.isot == "2010-10-15T23:00:11.000"


def test_date(mdi):
    assert mdi.date.isot == "2010-10-15T23:00:11.000"


def test_private_date_setters(mdi):
    _test_private_date_setters(mdi)


def test_instrument(mdi):
    """Tests the instrument property of the MDIMap object."""
    assert mdi.instrument == "MDI"


def test_waveunit(mdi):
    assert mdi.waveunit == "Angstrom"


def test_observer(mdi):
    assert mdi.observer_coordinate.frame.name == 'heliographic_stonyhurst'
    assert u.allclose(mdi.observer_coordinate.lat, Angle(mdi.meta['CRLT_OBS']*u.degree))
    assert u.allclose(mdi.observer_coordinate.radius, mdi.meta['DSUN_OBS']*u.m)


def test_carrington(mdi):
    assert u.allclose(mdi.carrington_longitude, Angle(mdi.meta['CRLN_OBS']*u.deg))
    assert u.allclose(mdi.carrington_latitude, Angle(mdi.meta['CRLT_OBS']*u.deg))


def test_unit(mdi):
    assert mdi.unit == u.dimensionless_unscaled


def test_synoptic_source(mdi_synoptic):
    assert isinstance(mdi_synoptic, MDISynopticMap)
    # Check that the WCS is valid
    with pytest.warns(SunpyMetadataWarning, match='Missing metadata for observer'):
        mdi_synoptic.wcs


def test_wcs(mdi, mdi_synoptic):
    # Smoke test that WCS is valid and can transform from pixels to world coordinates
    mdi.pixel_to_world(0*u.pix, 0*u.pix)
    with pytest.warns(SunpyMetadataWarning, match='Missing metadata for observer'):
        mdi_synoptic.pixel_to_world(0*u.pix, 0*u.pix)


def test_unit_synoptic(mdi_synoptic):
    assert mdi_synoptic.unit == u.G
    assert mdi_synoptic.unit == u.Unit("Mx/cm^2")
    assert mdi_synoptic.unit.to_string() == 'Mx / cm2'


def test_private_date_setters_synoptic(mdi_synoptic):
    _test_private_date_setters(mdi_synoptic)
