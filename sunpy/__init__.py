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
dev
    Experimental functions not meant for use in production
map
    Methods relating to the solar maps

"""
__version__ = 0.01

import sunpy.map
from sunpy.map import Map
from sunpy.map.MapCube import MapCube
from sunpy.map.CompositeMap import CompositeMap
from sunpy.dev import *
