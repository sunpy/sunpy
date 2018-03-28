"""
SunPy's LightCurve module provides a datatype for 1D time series data.

The objects also include data downloaders for their specific instruments, they
also support instantiation from files such as csv.
"""
from __future__ import absolute_import

__author__ = "Keith Hughitt"
__email__ = "keith.hughitt@nasa.gov"

import warnings
from sunpy.util.exceptions import SunpyDeprecationWarning as _SunpyDeprecationWarning

warnings.warn("As of v0.8.0, the `sunpy.lightcurve` module is deprecated and will be "
              "removed in a future version. Use `sunpy.timeseries` or "
              "`sunpy.map` for coordinate transformations.",
              _SunpyDeprecationWarning)

from sunpy.lightcurve.lightcurve import LightCurve
from sunpy.lightcurve.sources.eve import *
from sunpy.lightcurve.sources.goes import *
from sunpy.lightcurve.sources.noaa import *
from sunpy.lightcurve.sources.lyra import *
from sunpy.lightcurve.sources.logical import *
from sunpy.lightcurve.sources.norh import *
from sunpy.lightcurve.sources.rhessi import *
from sunpy.lightcurve.sources.fermi_gbm import *
