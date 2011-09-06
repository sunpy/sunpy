from __future__ import absolute_import

"""
SunPy
=====

An open-source Python library for Solar Physics data analysis.

Classes
-------
Sun
    Solar-related constants
Map
    A spatially-aware 2d data array
MapCube
    A spatially-aware 3d data array


Available subpackages
---------------------
map
    Methods relating to the solar maps

"""
#pylint: disable=W0401,W0622

__version__ = 0.01

import sunpy.map
import sunpy.data.sample
from sunpy.map import Map
from sunpy.map.mapcube import MapCube
from sunpy.map.compositemap import CompositeMap
from sunpy.cm import *

# Sample data
from sunpy.data.sample import AIA_171_IMAGE, RHESSI_IMAGE
