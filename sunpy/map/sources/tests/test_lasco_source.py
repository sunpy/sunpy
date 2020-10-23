"""Test cases for SOHO Map subclasses.
This particular test file pertains to LASCOMap.
@Author: Pritish C. (VaticanCameos)
"""

import os
import glob

import numpy as np
import pytest

import sunpy.data.test
from sunpy.map import Map
from sunpy.map.sources.soho import LASCOMap
from sunpy.tests.helpers import skip_glymur

path = sunpy.data.test.rootdir


@pytest.fixture
def lasco():
    fitspath = glob.glob(os.path.join(path, "lasco_c2_25299383_s.fts"))
    return Map(fitspath)


@pytest.fixture
def lasco_helioviewer():
    jp2path = glob.glob(os.path.join(
        path, "2013_05_13__16_54_06_137__SOHO_LASCO_C3_white-light.jp2"))
    return Map(jp2path)


# LASCO Tests


def test_fitstoLASCO(lasco):
    """Tests the creation of LASCOMap using FITS."""
    assert isinstance(lasco, LASCOMap)


def test_is_datasource_for(lasco):
    """Test the is_datasource_for method of LASCOMap.
    Note that header data to be provided as an argument
    can be a MetaDict object."""
    assert lasco.is_datasource_for(lasco.data, lasco.meta)


def test_measurement(lasco):
    """Tests the measurement property of the LASCOMap object."""
    assert lasco.measurement == "white-light"


def test_observatory(lasco):
    """Tests the observatory property of the LASCOMap object."""
    assert lasco.observatory == "SOHO"


def test_norm_clip(lasco):
    # Tests that the default normalizer has clipping disabled
    assert not lasco.plot_settings['norm'].clip


@skip_glymur
def test_helioviewer_rotation(lasco, lasco_helioviewer):
    """Tests that rotation metadata is correctly removed
    for JPEG2000 images provided by Helioviewer.org."""
    np.testing.assert_allclose(lasco.rotation_matrix,
                               [[0.999966, -0.008296], [0.008296, 0.999966]], rtol=1e-6)
    np.testing.assert_array_equal(lasco_helioviewer.rotation_matrix, [[1., 0.], [0., 1.]])
