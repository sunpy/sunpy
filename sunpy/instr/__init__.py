import warnings

from sunpy.util.exceptions import SunpyDeprecationWarning

warnings.warn("sunpy.instr is deprecated and will be removed in sunpy 3.1. "
              "Please use the affiliated sunkit-instruments package.",
              SunpyDeprecationWarning)
