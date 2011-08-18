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
#pylint: disable=W0401

__version__ = 0.01

import map
import data.sample
from map import Map
from map.mapcube import MapCube
from map.compositemap import CompositeMap
from dev import *
from cm import *

# Sample data
from data.sample import AIA_171_IMAGE, RHESSI_HSI
