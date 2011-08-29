"""SunPy sample data files"""
__author__ = "Keith Hughitt"
__email__ = "keith.hughitt@nasa.gov"

import sunpy
import os

rootdir = os.path.join(os.path.dirname(sunpy.__file__), "data/sample") 

#
# AIA20110319_105400_0171.fits
#
AIA_171_IMAGE = os.path.join(rootdir, "AIA20110319_105400_0171.fits")

#
# hsi_image_20110110_072142_072342_1__4_10_10_15_15_30keV_clean.fits
#
RHESSI_HSI = os.path.join(rootdir, "hsi_image_20110110_072142_072342_1__4_10_10_15_15_30keV_clean.fits")
