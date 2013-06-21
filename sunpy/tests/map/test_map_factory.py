# -*- coding: utf-8 -*-
"""
Created on Fri Jun 21 15:05:09 2013

@author: stuart
"""

import glob
import sunpy
import sunpy.map

filepath = "/home/stuart/SyncBox/Programming/SunPy/refactor/map_tests/"
a_list_of_many = glob.glob(filepath + "EIT/*")
a_fname = a_list_of_many[0]
#==============================================================================
# Map Factory Tests
#==============================================================================
class TestMap:
    #Test making a GenericMap
    
    def test_mapcube(self):
        #Test making a MapCube
        cube = sunpy.Map(*a_list_of_many, cube=True)
        assert isinstance(cube, sunpy.map.MapCube)
    
    def test_composite(self):
        #Test making a CompositeMap
        pass
    
    def test_patterns(self):
        ## Test different Map pattern matching ##
        # Data-header pair in a tuple
        
        # Data-header pair not in a tuple
         
        # File name
        eitmap = sunpy.Map(a_fname)
        assert isinstance(eitmap, sunpy.map.GenericMap)
        # Directory
        maps = sunpy.Map(filepath + "EIT/")
        assert isinstance(maps, list)
        assert ([isinstance(amap,sunpy.map.GenericMap) for amap in maps])
        # Glob
        maps = sunpy.Map(filepath + "EIT/*")
        assert isinstance(maps, list)
        assert ([isinstance(amap,sunpy.map.GenericMap) for amap in maps])
        # Already a Map
        amap = sunpy.Map(maps[0])
        assert isinstance(amap, sunpy.map.GenericMap)
        # A URL
        
        # A list of filenames
        #maps = sunpy.Map(a_list_of_many)
        #assert isinstance(maps, list)
        #assert ([isinstance(amap,sunpy.map.GenericMap) for amap in maps])
    
    def test_save(self):
        #Test save out
        pass
    
#==============================================================================
# Sources Tests
#==============================================================================
    def test_sdo(self):
        #Test an AIAMap & HMIMap
        aia = sunpy.Map(sunpy.AIA_171_IMAGE)
        assert isinstance(aia,sunpy.map.sources.AIAMap)
        
    def test_soho(self):
        #Test EITMap, LASCOMap & MDIMap
        eit = sunpy.Map(filepath + "EIT/efz20040301.000010.fits")
        assert isinstance(eit,sunpy.map.sources.EITMap)
    
        lasco = sunpy.Map(filepath + "lasco_c2_25299383.fts")
        assert isinstance(lasco,sunpy.map.sources.LASCOMap)
        
        mdi_c = sunpy.Map(filepath + "mdi_fd_Ic_6h_01d.5871.0000.fits")
        assert isinstance(mdi_c,sunpy.map.sources.MDIMap)
        
        mdi_m = sunpy.Map(filepath + "mdi_fd_M_96m_01d.5874.0005.fits")
        assert isinstance(mdi_m,sunpy.map.sources.MDIMap)
        
    def test_stereo(self):    
        #Test EUVIMap & CORMap
        euvi = sunpy.Map(filepath + "euvi_20090615_000900_n4euA.fts")
        assert isinstance(euvi,sunpy.map.sources.EUVIMap)
        
        cor = sunpy.Map(filepath + "cor1_20090615_000500_s4c1A.fts")
        assert isinstance(cor,sunpy.map.sources.CORMap)
        
    def test_rhessi(self):    
        #Test RHESSIMap
        rhessi = sunpy.Map(sunpy.RHESSI_IMAGE)
        assert isinstance(rhessi,sunpy.map.sources.RHESSIMap)
        
        #Test SWAPMap
        
        #Test XRTMap
        
        #Test SXTMap