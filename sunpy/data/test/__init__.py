"""SunPy test data files"""
from __future__ import absolute_import
import sunpy
import os

__author__ = "Keith Hughitt"
__email__ = "keith.hughitt@nasa.gov"

rootdir = os.path.join(os.path.dirname(sunpy.__file__), "data", "test") 

#
# EVE
#
EVE_LEVEL0_CSV = os.path.join(rootdir, "LATEST_EVE_L0CS_DIODES_1m.txt")
EVE_AVERAGES_CSV = os.path.join(rootdir, "EVE_He_II_304_averages.csv")
