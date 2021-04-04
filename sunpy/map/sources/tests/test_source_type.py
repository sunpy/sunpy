import os
import glob

import pytest

from astropy.visualization import LinearStretch

import sunpy.data.test
from sunpy.map import Map
from sunpy.map.sources.source_type import from_helioviewer_project, source_stretch
from sunpy.tests.helpers import skip_glymur

path = sunpy.data.test.rootdir


@pytest.fixture
def hvaia():
    fitspath = glob.glob(os.path.join(path, "aia_171_level1.fits"))
    return Map(fitspath)


@pytest.fixture
def hvjp2():
    jp2path = glob.glob(os.path.join(path, "2013_06_24__17_31_30_84__SDO_AIA_AIA_193.jp2"))
    return Map(jp2path)


@skip_glymur
def test_from_helioviewer_project(hvaia, hvjp2):
    """
    Tests if we are able to determine if a file is from the Helioviewer Project
    or not.
    """
    assert not from_helioviewer_project(hvaia.meta)
    assert from_helioviewer_project(hvjp2.meta)


@skip_glymur
def test_source_stretch(hvaia, hvjp2):
    """
    Tests that the correct stretch function is returned.
    """
    aia_fits_stretch = hvaia.plot_settings['norm'].stretch
    assert source_stretch(hvaia.meta, aia_fits_stretch) is aia_fits_stretch
    assert isinstance(source_stretch(hvjp2.meta, aia_fits_stretch), LinearStretch)
