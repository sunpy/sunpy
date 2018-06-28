"""Datasource-specific classes

This is where datasource specific logic is implemented. Each mission should
have its own file with one or more classes defined. Typically, these classes
will be subclasses of the :mod`sunpy.map.Map` class.
"""
from __future__ import division, print_function, absolute_import

from .sdo import AIAMap, HMIMap
from .iris import SJIMap
from .soho import EITMap, MDIMap, LASCOMap
from .trace import TRACEMap
from .hinode import SOTMap, XRTMap
from .proba2 import SWAPMap
from .rhessi import RHESSIMap
from .stereo import HIMap, CORMap, EUVIMap
from .yohkoh import SXTMap
from .source_type import source_stretch, from_helioviewer_project
from ..map_factory import Map

__all__ = ['XRTMap', 'SOTMap', 'SWAPMap', 'RHESSIMap', 'AIAMap', 'HMIMap',
           'EITMap', 'LASCOMap', 'MDIMap', 'EUVIMap', 'CORMap', 'HIMap',
           'SXTMap', 'SJIMap', 'TRACEMap', 'source_stretch']
