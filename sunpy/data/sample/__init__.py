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
# hsi_image_20101016_191218.fits
#
RHESSI_IMAGE = os.path.join(rootdir, "hsi_image_20101016_191218.fits")
