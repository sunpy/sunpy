"""Datasource-specific classes

This is where datasource specific logic is implemented. Each mission should
have its own file with one or more classes defined. Typically, these classes
will be subclasses of the :mod`sunpy.map.Map` class.
"""
from __future__ import absolute_import, division, print_function
__all__ = ['XRTMap', 'SOTMap', 'SWAPMap', 'RHESSIMap', 'AIAMap', 'HMIMap',
           'EITMap', 'LASCOMap', 'MDIMap', 'EUVIMap', 'CORMap', 'HIMap',
           'SXTMap', 'SJIMap', 'TRACEMap']

from .. map_factory import Map

from .hinode import XRTMap, SOTMap

from .proba2 import SWAPMap

from .rhessi import RHESSIMap

from .sdo import AIAMap, HMIMap

from .soho import EITMap, LASCOMap, MDIMap

from .stereo import EUVIMap, CORMap, HIMap

from .yohkoh import SXTMap

from .iris import SJIMap

from .trace import TRACEMap
