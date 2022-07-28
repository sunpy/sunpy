import numpy as np
import pytest

import astropy.units as u
from astropy.coordinates import Angle

from sunpy.data.test import get_dummy_map_from_header, get_test_filepath, test_data_filenames
from sunpy.map import Map
from sunpy.map.sources.soho import EITMap, LASCOMap, MDIMap, MDISynopticMap
from sunpy.tests.helpers import skip_glymur
from sunpy.time import parse_time
from sunpy.util.exceptions import SunpyMetadataWarning

efz_header_list = [f for f in test_data_filenames() if 'efz' in f.name and '.header' in f.name]
lasco_header_list = ["lasco_c2_25299383_s.header", "lasco_c3.header", ]


@pytest.fixture(scope="module", params=efz_header_list)
def eit_map(request):
    return get_dummy_map_from_header(request.param)


@pytest.fixture(scope="module", params=lasco_header_list, ids=['C2', 'C3'])
def lasco_map(request):
    return get_dummy_map_from_header(get_test_filepath(request.param))


@pytest.fixture
def lasco_helioviewer():
    return Map(get_test_filepath("2013_05_13__16_54_06_137__SOHO_LASCO_C3_white-light.jp2"))


@pytest.fixture
def mdi():
    return get_dummy_map_from_header(get_test_filepath("mdi.fd_Ic.20101015_230100_TAI.data.header"))


@pytest.fixture
def mdi_synoptic():
    return get_dummy_map_from_header(get_test_filepath('mdi_synoptic.header'))


def test_eit_map(eit_map):
    assert isinstance(eit_map, EITMap)


def test_is_datasource_for_eit(eit_map):
    assert eit_map.is_datasource_for(eit_map.data, eit_map.meta)


def test_observatory_eit(eit_map):
    assert eit_map.observatory == "SOHO"


def test_measurement_eit(eit_map):
    assert eit_map.measurement.value in [195, 171]


def test_rsun_eit(eit_map):
    assert u.allclose(eit_map.rsun_obs, 979.0701*u.arcsec)


def test_norm_clip_eit(eit_map):
    # Tests that the default normalizer has clipping disabled
    assert not eit_map.plot_settings['norm'].clip


def test_wcs_eit(eit_map):
    # Smoke test that WCS is valid and can transform from pixels to world coordinates
    eit_map.pixel_to_world(0*u.pix, 0*u.pix)


def test_lasco_map(lasco_map):
    assert isinstance(lasco_map, LASCOMap)


def test_is_datasource_for_lasco(lasco_map):
    assert lasco_map.is_datasource_for(lasco_map.data, lasco_map.meta)


def test_measurement_lasco(lasco_map):
    assert lasco_map.measurement == "white-light"


def test_wavelength_lasco(lasco_map):
    assert lasco_map.wavelength is None


def test_date_lasco(lasco_map):
    assert lasco_map.date == parse_time(
        {'C2': '2009-02-28T00:05:33.380',
         'C3': '2002-05-21T00:18:06.516'}[lasco_map.detector])


def test_nickname_lasco(lasco_map):
    assert lasco_map.nickname == {'C2': 'LASCO-C2 Orange',
                                  'C3': 'LASCO-C3 Clear'}[lasco_map.detector]


def test_observatory_lasco(lasco_map):
    assert lasco_map.observatory == "SOHO"


def test_norm_clip(lasco_map):
    # Tests that the default normalizer has clipping disabled
    assert not lasco_map.plot_settings['norm'].clip


@skip_glymur
def test_helioviewer_rotation_lasco(lasco_map, lasco_helioviewer):
    """
    Tests that rotation metadata is correctly removed for JPEG2000 images
    provided by Helioviewer.org.
    """
    rmatrix = {'C2': [[0.999966, -0.008296], [0.008296, 0.999966]],
               'C3': [[1, 0], [0, 1]]}[lasco_map.detector]
    np.testing.assert_allclose(lasco_map.rotation_matrix, rmatrix, rtol=1e-6)
    np.testing.assert_array_equal(lasco_helioviewer.rotation_matrix, [[1., 0.], [0., 1.]])


def test_wcs_lasco(lasco_map):
    # Smoke test that WCS is valid and can transform from pixels to world coordinates
    with pytest.warns(SunpyMetadataWarning, match='Missing metadata for observer'):
        lasco_map.pixel_to_world(0*u.pix, 0*u.pix)


def test_mdi_map(mdi):
    assert isinstance(mdi, MDIMap)


def test_is_datasource_for_mdi(mdi):
    assert mdi.is_datasource_for(mdi.data, mdi.meta)


def test_observatory_mdi(mdi):
    assert mdi.observatory == "SOHO"


def test_instrument_mdi(mdi):
    assert mdi.instrument == "MDI"


def test_waveunit_mdi(mdi):
    assert mdi.waveunit == "Angstrom"


def test_observer_mdi(mdi):
    assert mdi.observer_coordinate.frame.name == 'heliographic_stonyhurst'
    assert u.allclose(mdi.observer_coordinate.lat, Angle(mdi.meta['CRLT_OBS']*u.degree))
    assert u.allclose(mdi.observer_coordinate.radius, mdi.meta['DSUN_OBS']*u.m)


def test_carrington_mdi(mdi):
    assert u.allclose(mdi.carrington_longitude, Angle(mdi.meta['CRLN_OBS']*u.deg))
    assert u.allclose(mdi.carrington_latitude, Angle(mdi.meta['CRLT_OBS']*u.deg))


def test_unit_mdi(mdi):
    assert mdi.unit == u.dimensionless_unscaled


def test_synoptic_source(mdi_synoptic):
    assert isinstance(mdi_synoptic, MDISynopticMap)
    # Check that the WCS is valid
    with pytest.warns(SunpyMetadataWarning, match='Missing metadata for observer'):
        mdi_synoptic.wcs


def test_wcs_mdi(mdi, mdi_synoptic):
    # Smoke test that WCS is valid and can transform from pixels to world coordinates
    mdi.pixel_to_world(0*u.pix, 0*u.pix)
    with pytest.warns(SunpyMetadataWarning, match='Missing metadata for observer'):
        mdi_synoptic.pixel_to_world(0*u.pix, 0*u.pix)


def test_unit_synoptic(mdi_synoptic):
    assert mdi_synoptic.unit == u.G
    assert mdi_synoptic.unit == u.Unit("Mx/cm^2")
    assert mdi_synoptic.unit.to_string() == 'Mx / cm2'
