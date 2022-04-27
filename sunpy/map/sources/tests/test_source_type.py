"""
Test cases distinguishing the source types.
"""
import pytest

from astropy.visualization import LinearStretch

from sunpy.data.test import get_test_filepath
from sunpy.map import Map
from sunpy.map.sources.source_type import from_helioviewer_project, source_stretch
from sunpy.tests.helpers import skip_glymur


@pytest.fixture
def hvjp2():
    return Map(get_test_filepath("2013_06_24__17_31_30_84__SDO_AIA_AIA_193.jp2"))


@skip_glymur
def test_from_helioviewer_project(aia171_test_map, hvjp2):
    """Tests if we are able to determine if a file is from the Helioviewer
    Project or not."""
    assert not from_helioviewer_project(aia171_test_map.meta)
    assert from_helioviewer_project(hvjp2.meta)


@skip_glymur
def test_source_stretch(aia171_test_map, hvjp2):
    """
    Tests that the correct stretch function is returned.
    """
    aia_fits_stretch = aia171_test_map.plot_settings['norm'].stretch
    assert source_stretch(aia171_test_map.meta, aia_fits_stretch) is aia_fits_stretch
    assert isinstance(source_stretch(hvjp2.meta, aia_fits_stretch), LinearStretch)
