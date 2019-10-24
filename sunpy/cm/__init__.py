import sys as _sys

# Deprecate the module without polluting the namespace
import warnings as _w
from sunpy.util.exceptions import SunpyDeprecationWarning as _SDW

_w.warn("The sunpy.cm module has been moved in 1.1 to sunpy.visualization.colormaps, it will be removed in 2.1",
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
