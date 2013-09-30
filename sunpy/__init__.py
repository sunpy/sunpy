"""
SunPy
=====

An open-source Python library for Solar Physics data analysis.

Web Links
---------
Homepage: http://www.sunpy.org
Documentation: http://sunpy.readthedocs.org/en/latest/index.html

Available subpackages
---------------------
cm
    Solar Physics specific color maps.
image
    Generic methods to work with image data or maps.
instr
    subpackages spcific to missions or instruments.
lightcurve
    subpackage for working with 1D data sets like lightcurves.
map
    subpackage for working with 2D and 3D data sets of images or sequences of 
    images.
net
    Routines for obtaining data from the internet.
spectra
    subpackage for working with 2D spectra datatypes
sun
    Contains astronomical and physical constants.
time
    Contains time related constants and methods.
wcs
     The WCS package provides functions to parse a World Coordinate System(WCS)
coordinates for solar images as well as convert between various solar
coordinate systems.

"""
from __future__ import absolute_import

__version__ = 0.3

import sunpy.map
import sunpy.sun
import sunpy.lightcurve

from sunpy.map.map_factory import Map, make_map

from sunpy.util.config import load_config, print_config
from sunpy.cm import *

# Sample data
from sunpy.data.sample import (AIA_171_IMAGE, RHESSI_IMAGE, EIT_195_IMAGE, 
                               RHESSI_EVENT_LIST, CALLISTO_IMAGE)

# Load user configuration
config = load_config()
