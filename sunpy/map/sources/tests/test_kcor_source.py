import os
import glob

import pytest

import astropy.units as u

import sunpy.data.test
from sunpy.map import Map
from sunpy.map.sources.mlso import KCorMap


@pytest.fixture()
def kcor():
    """
    Creates an KCorMap from a FITS file.
    """
    path = sunpy.data.test.rootdir
    fitspath = glob.glob(os.path.join(path, "20181209_180305_kcor_l1.5_rebinned.fits"))
    return Map(fitspath)


def test_kcormap_creation(kcor):
    """
    Tests the creation of KCorMap using FITS.
    """
    assert isinstance(kcor, KCorMap)


def test_is_datasource_for(kcor):
    """
    Test the is_datasource_for method of KCorMap.
    """
    assert kcor.is_datasource_for(kcor.data, kcor.meta)


def test_measurement(kcor):
    """
    Tests the measurement property of KCorMap.
    """
    assert kcor.measurement == 735 * u.nm


def test_observatory(kcor):
    """
    Tests the observatory property of KCorMap.
    """
    assert kcor.observatory == "MLSO"


def test_norm_clip(kcor):
    """
    Tests that the default normalizer has clipping disabled.
    """
    assert not kcor.plot_settings['norm'].clip
