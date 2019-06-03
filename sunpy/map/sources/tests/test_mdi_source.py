"""Test cases for SOHO Map subclasses.
This particular test file pertains to MDIMap.
@Author: Pritish C. (VaticanCameos)
"""

import os
import glob

import pytest

import astropy.units as u

from sunpy.coordinates import frames
from sunpy.map.sources.soho import MDIMap
from sunpy.map import Map
import sunpy.data.test


@pytest.fixture
def mdi():
    path = sunpy.data.test.rootdir
    fitspath = glob.glob(os.path.join(path, "mdi_fd_Ic_6h_01d.5871.0000_s.fits"))
    return Map(fitspath)


# MDI Tests
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


def test_measurement(mdi):
    """Tests the measurement property of the MDIMap object."""
    assert mdi.measurement == "continuum"


def test_waveunit(mdi):
    assert mdi.waveunit == "Angstrom"


def test_observer(mdi):
    assert isinstance(mdi.observer_coordinate.frame, frames.HeliographicStonyhurst)
    assert u.allclose(mdi.observer_coordinate.lat, -5.774028172878*u.deg)
    assert u.allclose(mdi.observer_coordinate.lon, -0.10181577*u.deg)
    assert u.allclose(mdi.observer_coordinate.radius, 0.9739569156244*u.AU)
