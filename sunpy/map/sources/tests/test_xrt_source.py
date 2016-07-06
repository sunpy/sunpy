"""Test cases for HINODE Map subclasses.
This particular test file pertains to XRTMap.
@Author: Pritish C. (VaticanCameos)
"""

import os
import glob
import numpy as np
from matplotlib import colors

import astropy.units as u

from sunpy.map.sources.hinode import XRTMap
from sunpy.map import Map
import sunpy.data.test

path = sunpy.data.test.rootdir
fitspath = glob.glob(os.path.join(path, "HinodeXRT.fits"))
xrt = Map(fitspath)


# XRT Tests
def test_fitstoXRT():
    """Tests the creation of XRTMap using FITS."""
    assert isinstance(xrt, XRTMap)


def test_is_datasource_for():
    """Test the is_datasource_for method of XRTMap.
    Note that header data to be provided as an argument
    can be a MetaDict object."""
    assert xrt.is_datasource_for(xrt.data, xrt.meta)


def test_observatory():
    """Tests the observatory property of the XRTMap object."""
    assert xrt.observatory == "Hinode"


def test_measurement():
    """Tests the measurement property of the XRTMap object."""
    measurement = xrt.filter_wheel1_measurements[5].replace("_", " ")
    measurement += '-' + xrt.filter_wheel2_measurements[1].replace("_", " ")
    assert xrt.measurement == measurement


def test_wheel_measurements():
    """Tests the filter_wheel_measurements objects present
    in the XRTMap object."""
    assert (xrt.filter_wheel1_measurements ==
            ["Al_med", "Al_poly", "Be_med", "Be_thin", "C_poly", "Open"])
    assert (xrt.filter_wheel2_measurements ==
            ["Open", "Al_mesh", "Al_thick", "Be_thick", "Gband", "Ti_poly"])
