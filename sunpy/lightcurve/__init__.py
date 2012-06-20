"""SunPy LightCurves"""
from __future__ import absolute_import

__author__ = "Keith Hughitt"
__email__ = "keith.hughitt@nasa.gov"

from sunpy.lightcurve.lightcurve import LightCurve
from sunpy.lightcurve.sources.sdo import *

def make_lightcurve(*args, **kwargs):
    """Processes one or more inputs and returns a LightCurve instance.
    
    Parameters
    ----------
    args : filepath
        The data source used to create the lightcurve object.
        
    Returns
    -------
    out : LightCurve
        Returns a LightCurve subclass instance
        
    Examples
    --------
    >>> import sunpy
    >>> sunpy.make_lightcurve("file.fts")
    """
    if len(args) is 0:
        raise TypeError("Invalid input.")

    return LightCurve.read(args[0])
