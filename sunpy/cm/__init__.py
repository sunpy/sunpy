import sys as _sys

# Deprecate the module without polluting the namespace
import warnings as _w
from sunpy.util.exceptions import SunpyDeprecationWarning as _SDW

_w.warn("The functionality of the sunpy.cm module is now in sunpy.visualization.colormaps as of SunPy 1.1. "
        "The ability to import sunpy.cm as a module will be removed in a future version of SunPy.",
        category=_SDW)

# Import docstring from the new module
import sunpy.visualization.colormaps

__doc__ = sunpy.visualization.colormaps.__doc__

# Trick Python into being able to do from sunpy.cm.cm import sdoaia171
# This does mean these two are imported by default where only cm was before.
from sunpy.visualization.colormaps import color_tables
from sunpy.visualization.colormaps import cm

_sys.modules['sunpy.cm.cm'] = sunpy.visualization.colormaps.cm
_sys.modules['sunpy.cm.color_tables'] = sunpy.visualization.colormaps.color_tables

# Import the two original functions in this namespace and use __all__ to keep
# import * behaviour the same.
from sunpy.visualization.colormaps.cm import *

__all__ = ['show_colormaps', 'cmlist']
