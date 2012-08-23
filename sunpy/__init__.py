"""
SunPy
=====

An open-source Python library for Solar Physics data analysis.

Classes
-------
Map
    A spatially-aware 2d data array
MapCube
    A spatially-aware 3d data array


Available subpackages
---------------------
map
    Methods relating to the solar maps

"""
from __future__ import absolute_import

__version__ = 0.1

import sunpy.map
import sunpy.sun
import sunpy.lightcurve
from sunpy.map import make_map
from sunpy.map import read_header
from sunpy.map.map import Map
from sunpy.map.header import MapHeader
from sunpy.map.mapcube import MapCube
from sunpy.map.compositemap import CompositeMap
from sunpy.util.config import load_config, print_config
from sunpy.cm import *

# Sample data
from sunpy.data.sample import (AIA_171_IMAGE, RHESSI_IMAGE, EIT_195_IMAGE, 
                               RHESSI_EVENT_LIST, CALLISTO_IMAGE)

# Load user configuration
config = load_config()

# Keep track of created data types
from sunpy.util.data_manager import DataManager
_data_manager = DataManager()
