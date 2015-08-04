"""Datasource-specific classes

This is where datasource specific logic is implemented. Each mission should
have its own file with one or more classes defined. Typically, these classes
will be subclasses of the :mod`sunpy.map.Map` class.
"""
__all__ = ['XRTMap', 'SOTMap', 'SWAPMap', 'RHESSIMap', 'AIAMap', 'HMIMap',
           'EITMap', 'LASCOMap', 'MDIMap', 'EUVIMap', 'CORMap', 'HIMap',
           'SXTMap', 'SJIMap', 'TRACEMap']

from .. map_factory import Map

from hinode import XRTMap, SOTMap
Map.register(XRTMap, XRTMap.is_datasource_for)
Map.register(SOTMap, SOTMap.is_datasource_for)

from proba2 import SWAPMap
Map.register(SWAPMap, SWAPMap.is_datasource_for)

from rhessi import RHESSIMap
Map.register(RHESSIMap, RHESSIMap.is_datasource_for)

from sdo import AIAMap, HMIMap
Map.register(AIAMap, AIAMap.is_datasource_for)
Map.register(HMIMap, HMIMap.is_datasource_for)

from soho import EITMap, LASCOMap, MDIMap
Map.register(EITMap, EITMap.is_datasource_for)
Map.register(LASCOMap, LASCOMap.is_datasource_for)
Map.register(MDIMap, MDIMap.is_datasource_for)

from stereo import EUVIMap, CORMap, HIMap
Map.register(EUVIMap, EUVIMap.is_datasource_for)
Map.register(CORMap, CORMap.is_datasource_for)
Map.register(HIMap, HIMap.is_datasource_for)

from yohkoh import SXTMap
Map.register(SXTMap, SXTMap.is_datasource_for)

from iris import SJIMap
Map.register(SJIMap, SJIMap.is_datasource_for)

from trace import TRACEMap
Map.register(TRACEMap, TRACEMap.is_datasource_for)
