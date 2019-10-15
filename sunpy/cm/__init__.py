from sunpy.util.decorators import deprecated as _d

_d("1.1",
   "The sunpy.cm module has been moved.",
   alternative="sunpy.visualization.colormaps")

from sunpy.visualization.colormaps import cm
from sunpy.visualization.colormaps import color_tables
from sunpy.visualization.colormaps.cm import *
