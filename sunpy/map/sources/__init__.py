"""Datasource-specific classes

This is where datasource specific logic is implemented. Each mission should
have its own file with one or more classes defined. Typically, these classes
will be subclasses of the :mod`sunpy.map.Map` class.
"""
__all__ = ['XRTMap', 'SOTMap', 'SWAPMap', 'RHESSIMap', 'AIAMap', 'HMIMap',
           'EITMap', 'LASCOMap', 'MDIMap', 'EUVIMap', 'CORMap', 'HIMap',
           'SXTMap', 'SJIMap', 'TRACEMap', 'KCorMap', 'SUVIMap',
           'source_stretch']

from ..map_factory import Map
from .hinode import SOTMap, XRTMap
from .iris import SJIMap
from .mlso import KCorMap
from .proba2 import SWAPMap
from .rhessi import RHESSIMap
from .sdo import AIAMap, HMIMap
from .soho import EITMap, LASCOMap, MDIMap
from .source_type import from_helioviewer_project, source_stretch
from .stereo import CORMap, EUVIMap, HIMap
from .suvi import SUVIMap
from .trace import TRACEMap
from .yohkoh import SXTMap
