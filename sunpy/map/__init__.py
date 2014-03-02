"""SunPy Maps"""
from __future__ import absolute_import

__author__ = "Keith Hughitt"
__email__ = "keith.hughitt@nasa.gov"


from sunpy.map.mapbase import GenericMap

from sunpy.map.header import MapMeta
from . mapcube import MapCube
from . compositemap import CompositeMap

from sunpy.map.map_factory import Map
#from sunpy.map.sources import *

from sunpy.instr.hinode import XRTMap
Map.register(XRTMap, XRTMap.is_datasource_for)

from sunpy.instr.iris import IRISMap
Map.register(IRISMap, IRISMap.is_datasource_for)

from sunpy.instr.proba2 import SWAPMap
Map.register(SWAPMap, SWAPMap.is_datasource_for)

from sunpy.instr.rhessi import RHESSIMap
Map.register(RHESSIMap, RHESSIMap.is_datasource_for)

from sunpy.instr.sdo import AIAMap, HMIMap
Map.register(AIAMap, AIAMap.is_datasource_for)
Map.register(HMIMap, HMIMap.is_datasource_for)

from sunpy.instr.soho import EITMap, LASCOMap, MDIMap
Map.register(EITMap, EITMap.is_datasource_for)
Map.register(LASCOMap, LASCOMap.is_datasource_for)
Map.register(MDIMap, MDIMap.is_datasource_for)

from sunpy.instr.stereo import EUVIMap, CORMap
Map.register(EUVIMap, EUVIMap.is_datasource_for)
Map.register(CORMap, CORMap.is_datasource_for)

from sunpy.instr.trace import TRACEMap
Map.register(TRACEMap, TRACEMap.is_datasource_for)

from sunpy.instr.yohkoh import SXTMap
Map.register(SXTMap, SXTMap.is_datasource_for)

