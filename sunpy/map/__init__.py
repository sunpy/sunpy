"""SunPy Maps"""
from __future__ import absolute_import

__author__ = "Keith Hughitt"
__email__ = "keith.hughitt@nasa.gov"


from sunpy.map.mapbase import GenericMap

from . mapcube import MapCube
from . mapsequence import MapSequence
from . compositemap import CompositeMap

from sunpy.map.map_factory import Map
from sunpy.map import sources
