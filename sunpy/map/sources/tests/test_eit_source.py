"""Test cases for SOHO Map subclasses.
This particular test file pertains to EITMap.
@Author: Pritish C. (VaticanCameos)
"""

import os
import glob

import pytest

import astropy.units as u

import sunpy.data.test
from sunpy.map import Map
from sunpy.map.sources.soho import EITMap

path = sunpy.data.test.rootdir
fitslist = glob.glob(os.path.join(path, "EIT", "*"))


@pytest.fixture(scope="module", params=fitslist)
def createEIT(request):
    """Creates an EITMap from a FITS file."""
    return Map(request.param)


# EIT Tests
def test_fitstoEIT(createEIT):
    """Tests the creation of EITMap using FITS."""
    assert isinstance(createEIT, EITMap)


def test_is_datasource_for(createEIT):
    """Test the is_datasource_for method of EITMap.
    Note that header data to be provided as an argument
    can be a MetaDict object."""
    assert createEIT.is_datasource_for(createEIT.data, createEIT.meta)


def test_observatory(createEIT):
    """Tests the observatory property of the EITMap object."""
    assert createEIT.observatory == "SOHO"


def test_measurement(createEIT):
    """Tests the measurement property of the EITMap object."""
    assert createEIT.measurement.value in [195, 171]


def test_rsun(createEIT):
    """Tests the measurement property of the EITMap object."""
    assert u.allclose(createEIT.rsun_obs, 979.0701*u.arcsec)


def test_norm_clip(createEIT):
    # Tests that the default normalizer has clipping disabled
    assert not createEIT.plot_settings['norm'].clip
