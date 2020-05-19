"""Test cases for HINODE Map subclasses.
This particular test file pertains to SOTMap.
@Author: Pritish C. (VaticanCameos)
"""

import os
import glob

import pytest

import astropy.units as u

import sunpy.data.test
from sunpy.map import Map
from sunpy.map.sources.hinode import SOTMap
from sunpy.util.exceptions import SunpyUserWarning


@pytest.fixture
def sot():
    path = sunpy.data.test.rootdir
    fitspath = glob.glob(os.path.join(path, "HinodeSOT.fits"))
    with pytest.warns(SunpyUserWarning, match='This file contains more than 2 dimensions'):
        return Map(fitspath)


# SOT Tests
def test_fitstoSOT(sot):
    """Tests the creation of SOTMap using FITS."""
    assert isinstance(sot, SOTMap)


def test_is_datasource_for(sot):
    """Test the is_datasource_for method of SOTMap.
    Note that header data to be provided as an argument
    can be a MetaDict object."""
    assert sot.is_datasource_for(sot.data, sot.meta)


def test_observatory(sot):
    """Tests the observatory property of the SOTMap object."""
    assert sot.observatory == "Hinode"


def test_measurement(sot):
    """Tests the measurement property of the SOTMap object."""
    assert sot.measurement == 0 * u.one


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
