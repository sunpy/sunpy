"""
The `sunpy.net.dataretriever` submodule is a framework for downloading data
from "simple" web sources such as HTTP or FTP servers. Although it could be
used for more complex services as well. Following the example of `sunpy.map`
and `sunpy.timeseries` this module provides a base class
`~sunpy.net.dataretriever.GenericClient` from which specific services can
subclass. All these subclasses are then registered with the `sunpy.net.Fido` factory
class, so do not need to be called individually.
"""

from .client import *
from .sources.eve import *
from .sources.fermi_gbm import *
from .sources.goes import *
from .sources.gong import *
from .sources.lyra import *
from .sources.noaa import *
from .sources.norh import *
from .sources.rhessi import *
