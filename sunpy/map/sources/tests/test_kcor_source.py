"""
Test cases for KCor Map subclass.
"""

import os
import glob

import pytest

import astropy.units as u

from sunpy.map.sources.mlso import KCorMap
from sunpy.map import Map
import sunpy.data.test

path = sunpy.data.test.rootdir
fitspath = glob.glob(os.path.join(path, "20181209_180305_kcor_l1.5_rebinned.fits"))


@pytest.fixture()
def kcor():

    """Creates an KCorMap from a FITS file."""
    return Map(fitspath)


# KCor Tests
def test_kcormap_creation(kcor):
    """Tests the creation of KCorMap using FITS."""
    assert isinstance(kcor, KCorMap)


def test_is_datasource_for(kcor):
    """Test the is_datasource_for method of KCorMap.
    Note that header data to be provided as an argument
    can be a MetaDict object."""
    assert kcor.is_datasource_for(kcor.data, kcor.meta)


def test_measurement(kcor):
    """Tests the measurement property of the KCorMap object."""
    assert kcor.measurement == 735 * u.nm


def test_observatory(kcor):
    """Tests the observatory property of the KCorMap object."""
    assert kcor.observatory == "MLSO"
