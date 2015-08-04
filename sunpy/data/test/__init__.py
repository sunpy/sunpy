"""SunPy test data files"""
from __future__ import absolute_import
import sunpy
import os
import glob

__author__ = "Keith Hughitt"
__email__ = "keith.hughitt@nasa.gov"

rootdir = os.path.join(os.path.dirname(sunpy.__file__), "data", "test")

#
# EVE
#
EVE_LEVEL0_CSV = os.path.join(rootdir, "LATEST_EVE_L0CS_DIODES_1m.txt")
EVE_AVERAGES_CSV = os.path.join(rootdir, "EVE_He_II_304_averages.csv")

#
# JPEG2000 sample
#
AIA_193_JP2 = os.path.join(rootdir,
                           "2013_06_24__17_31_30_84__SDO_AIA_AIA_193.jp2")

#
# aiaprep() test Map
#
aia_171_level1 = os.path.join(rootdir, "aia_171_level1.fits")

file_list = glob.glob(os.path.join(rootdir, '*.[!p]*'))
