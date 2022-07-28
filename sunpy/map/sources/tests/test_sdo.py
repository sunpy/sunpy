import pytest

import astropy.units as u
from astropy.visualization.mpl_normalize import ImageNormalize

from sunpy.data.test import get_dummy_map_from_header, get_test_filepath
from sunpy.map import Map
from sunpy.map.sources.sdo import AIAMap, HMIMap, HMISynopticMap
from sunpy.tests.helpers import SKIP_GLYMUR
from sunpy.util.exceptions import SunpyMetadataWarning

__author__ = "It's a me, Mario!"


params = [get_test_filepath("aia_171_level1.fits")]
if not SKIP_GLYMUR:
    params += [get_test_filepath("2013_06_24__17_31_30_84__SDO_AIA_AIA_193.jp2")]


@pytest.fixture(scope="module", params=params)
def aia_map(request):
    return Map(request.param)


@pytest.fixture
def hmi_map():
    return Map(get_test_filepath('resampled_hmi.fits'))


@pytest.fixture
def hmi_synoptic():
    return get_dummy_map_from_header(get_test_filepath('hmi_synoptic.header'))


def test_aia_map(aia_map):
    assert isinstance(aia_map, AIAMap)


def test_is_datasource_for_aia(aia_map):
    """
    Tests the is_datasource_for method of AIAMap.
    """
    assert aia_map.is_datasource_for(aia_map.data, aia_map.meta)


def test_observatory_aia(aia_map):
    assert aia_map.observatory == "SDO"


def test_measurement_aia(aia_map):
    # aiaimg has 171, jp2path has 193.
    assert aia_map.measurement.value in [171, 193]


def test_norm_clip_aia(aia_map):
    # Tests that the default normalizer has clipping disabled
    assert not aia_map.plot_settings['norm'].clip


def test_new_instance_preserves_plot_settings_aia(aia_map):
    # Tests that the _new_instance method preserves the plot_settings
    # of the old instance. This is done on the AIA source as the AIAMap
    # constructor explicitly sets the cmap and norm and we want to test
    # that _new_instance persists the old custom plot_settings
    aia_map.plot_settings['norm'] = ImageNormalize(vmin=0.1, vmax=42)
    aia_map.plot_settings['cmap'] = 'inferno'
    new_aia_map = aia_map._new_instance(aia_map.data,
                                        aia_map.meta,
                                        plot_settings=aia_map.plot_settings)
    assert new_aia_map.plot_settings['norm'].vmin == aia_map.plot_settings['norm'].vmin
    assert new_aia_map.plot_settings['norm'].vmax == aia_map.plot_settings['norm'].vmax
    assert new_aia_map.plot_settings['cmap'] == aia_map.plot_settings['cmap']
    # If no plot settings are explicitly passed, the plot_settings should fall back to those
    # in the constructor
    new_aia_map = aia_map._new_instance(aia_map.data, aia_map.meta)
    assert new_aia_map.plot_settings['norm'].vmin is None
    assert new_aia_map.plot_settings['norm'].vmax is None
    assert new_aia_map.plot_settings['cmap'] == new_aia_map._get_cmap_name()


def test_wcs_aia(aia_map):
    # Smoke test that WCS is valid and can transform from pixels to world coordinates
    aia_map.pixel_to_world(0*u.pix, 0*u.pix)


def test_hmi_map(hmi_map):
    assert isinstance(hmi_map, HMIMap)


def test_is_datasource_for_hmi(hmi_map):
    assert hmi_map.is_datasource_for(hmi_map.data, hmi_map.meta)


def test_observatory_hmi(hmi_map):
    assert hmi_map.observatory == "SDO"


def test_measurement_hmi(hmi_map):
    assert hmi_map.measurement == "continuum"


def test_wcs_hmi(hmi_map):
    # Smoke test that WCS is valid and can transform from pixels to world coordinates
    hmi_map.pixel_to_world(0*u.pix, 0*u.pix)


def test_hmi_synoptic_map(hmi_synoptic):
    assert isinstance(hmi_synoptic, HMISynopticMap)


def test_is_datasource_for_synoptic(hmi_synoptic):
    assert hmi_synoptic.is_datasource_for(hmi_synoptic.data, hmi_synoptic.meta)


def test_observatory_synoptic(hmi_synoptic):
    assert hmi_synoptic.observatory == "SDO"


def test_measurement_synoptic(hmi_synoptic):
    assert hmi_synoptic.measurement == "carrington"


def test_date_uses_date_obs_synoptic(hmi_synoptic):
    # Check that accessing the date doesn't raise a warning.
    hmi_synoptic.date
    hmi_synoptic.meta['date-obs'] = hmi_synoptic.meta.pop('t_obs')
    assert hmi_synoptic.date is not None


def test_unit_synoptic(hmi_synoptic):
    assert hmi_synoptic.unit == u.G
    assert hmi_synoptic.unit == u.Unit("Mx/cm^2")
    assert hmi_synoptic.unit.to_string() == 'Mx / cm2'
    hmi_synoptic.meta['bunit'] = 'm'
    assert hmi_synoptic.unit == u.m


def test_wcs_synoptic(hmi_synoptic):
    # Smoke test that WCS is valid and can transform from pixels to world coordinates
    with pytest.warns(SunpyMetadataWarning, match='Missing metadata for observer'):
        hmi_synoptic.pixel_to_world(0*u.pix, 0*u.pix)
