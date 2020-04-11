"""SunPy Maps"""
import sys

from sunpy.map.mapbase import GenericMap  # isort:skip

from sunpy.map import sources
from sunpy.map.header_helper import *
from sunpy.map.map_factory import Map
from sunpy.map.maputils import *

from .compositemap import CompositeMap
from .mapsequence import MapSequence


# If the documentation is being built using Sphinx, append output HTML to the docstrings for the
#   GenericMap and MapSequence `quicklook()` methods.
# This code detects a Sphinx build via the simplistic approach of checking if Sphinx has been
#   imported, and thus there is the potential for false positives.
if 'sphinx' in sys.modules:
    GenericMap._append_quicklook_example_to_docstring()
    MapSequence._append_quicklook_example_to_docstring()
