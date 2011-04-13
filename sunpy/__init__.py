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
import dev
import data
import Sun
from Map import Map
from MapCube import MapCube

