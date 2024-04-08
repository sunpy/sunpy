"""
SunPy Map

isort:skip_file
"""
# Check if user has installed the map extras
from sunpy.util.sysinfo import _warn_missing_deps
_warn_missing_deps('map')

from sunpy.map.mapbase import *

from sunpy.map import sources
from sunpy.map.header_helper import *
from sunpy.map.map_factory import Map
from sunpy.map.maputils import *
from .compositemap import CompositeMap
from .mapsequence import MapSequence
