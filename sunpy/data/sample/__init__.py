"""SunPy sample data files"""
from __future__ import absolute_import
import sunpy
import os

__author__ = "Keith Hughitt"
__email__ = "keith.hughitt@nasa.gov"

rootdir = os.path.join(os.path.dirname(sunpy.__file__), "data", "sample") 

#
# AIA20110319_105400_0171.fits
#
AIA_171_IMAGE = os.path.abspath(
    os.path.join(rootdir, "AIA20110319_105400_0171.fits")
)

#
# hsi_image_20101016_191218.fits
#
RHESSI_IMAGE = os.path.abspath(
    os.path.join(rootdir, "hsi_image_20101016_191218.fits")
)

#
# eit_l1_20020625_100011.fits
#
EIT_195_IMAGE = os.path.abspath(
    os.path.join(rootdir, "eit_l1_20020625_100011.fits")
)

CALLISTO_IMAGE = os.path.abspath(
    os.path.join(rootdir, "BIR_20110922_103000_01.fit")
)

# A stacked calibrated event list from RHESSI
# hsi_calib_ev_20020220_1106_20020220_1106_25_40.fits
#
RHESSI_EVENT_LIST = os.path.abspath(
    os.path.join(rootdir, "hsi_calib_ev_20020220_1106_20020220_1106_25_40.fits")
)

#
# swap_lv1_20120101_001607.fits
#
SWAP_LEVEL1_IMAGE = os.path.join(rootdir, "swap_lv1_20120101_001607.fits")
