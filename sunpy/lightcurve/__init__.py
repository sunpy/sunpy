"""
SunPy's LightCurve module provides a datatype for 1D time series data.

The objects also include data downloaders for their specific instruments, they
also support instantiation from files such as csv.
"""
from __future__ import absolute_import

__author__ = "Keith Hughitt"
__email__ = "keith.hughitt@nasa.gov"

from sunpy.lightcurve.lightcurve import LightCurve
from sunpy.lightcurve.sources.eve import *
from sunpy.lightcurve.sources.goes import *
from sunpy.lightcurve.sources.noaa import *
from sunpy.lightcurve.sources.lyra import *
from sunpy.lightcurve.sources.logical import *
from sunpy.lightcurve.sources.norh import *
from sunpy.lightcurve.sources.rhessi import *
from sunpy.lightcurve.sources.fermi_gbm import *

#def make_lightcurve(*args, **kwargs):
#    """Processes one or more inputs and returns a LightCurve instance.
#
#    Parameters
#    ----------
#    args : filepath, url, or start and end dates
#        The input for a LightCurve object should either be a filepath, a URL,
#        or a date range to be queried for the particular instrument.
#
#    Returns
#    -------
#    out : LightCurve
#        Returns a LightCurve subclass instance
#
#    Examples
#    --------
#    >>> import sunpy
#    >>> sunpy.make_lightcurve("file.fts")
#    """
#    if len(args) is 0:
#        raise TypeError("Invalid input.")
#
#    return LightCurve.read(args[0])
