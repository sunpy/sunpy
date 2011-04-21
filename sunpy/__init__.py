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
data
    Methods relating to the parsing of common solar data formats

"""
__version__ = 0.01

import sunpy.data
import sunpy.Sun
from sunpy.data.map import Map
from sunpy.data.MapCube import MapCube
from sunpy.data.CompositeMap import CompositeMap
from sunpy.dev import *
