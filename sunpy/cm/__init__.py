import warnings as _w
from sunpy.util.exceptions import SunpyDeprecationWarning

_w.warn("The sunpy.cm module has been moved in 1.1 to sunpy.visualization.colormaps, it will be removed in 2.1", category=SunpyDeprecationWarning)

import sunpy.visualization.colormaps
from sunpy.visualization.colormaps import cm
from sunpy.visualization.colormaps import color_tables
from sunpy.visualization.colormaps.cm import *


__doc__ = sunpy.visualization.colormaps.__doc__
