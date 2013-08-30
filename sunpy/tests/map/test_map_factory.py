# -*- coding: utf-8 -*-
"""
Created on Fri Jun 21 15:05:09 2013

@author: stuart
"""
import os
import glob
import sunpy
import sunpy.map
import sunpy.data.test

filepath = sunpy.data.test.rootdir
a_list_of_many = glob.glob(os.path.join(filepath,"EIT")+'/*')
a_fname = a_list_of_many[0]
#==============================================================================
# Map Factory Tests
#==============================================================================
class TestMap:
    def test_mapcube(self):
        #Test making a MapCube
        cube = sunpy.Map(a_list_of_many, cube=True)
        assert isinstance(cube, sunpy.map.MapCube)
    
    def test_composite(self):
        #Test making a CompositeMap
        comp = sunpy.Map(sunpy.AIA_171_IMAGE, sunpy.RHESSI_IMAGE,
                         composite=True)
        assert isinstance(comp, sunpy.map.CompositeMap)
    
    def test_patterns(self):
        ## Test different Map pattern matching ##
        # File name
        eitmap = sunpy.Map(a_fname)
        assert isinstance(eitmap, sunpy.map.GenericMap)
        # Directory
        maps = sunpy.Map(os.path.join(filepath,"EIT"))
        assert isinstance(maps, list)
        assert ([isinstance(amap,sunpy.map.GenericMap) for amap in maps])
        # Glob
        maps = sunpy.Map(os.path.join(filepath,"EIT")+'/*')
        assert isinstance(maps, list)
        assert ([isinstance(amap,sunpy.map.GenericMap) for amap in maps])
        # Already a Map
        amap = sunpy.Map(maps[0])
        assert isinstance(amap, sunpy.map.GenericMap)
        # A URL
        amap = sunpy.Map("https://raw.github.com/sunpy/sunpy/master/sunpy/data/sample/AIA20110319_105400_0171.fits")
        assert isinstance(amap, sunpy.map.GenericMap)
        # A list of filenames
        maps = sunpy.Map(a_list_of_many)
        assert isinstance(maps, list)
        assert ([isinstance(amap,sunpy.map.GenericMap) for amap in maps])
        # Data-header pair in a tuple
        pair_map = sunpy.Map((amap.data, amap.meta))
        assert isinstance(pair_map, sunpy.map.GenericMap)
        # Data-header pair not in a tuple
        pair_map = sunpy.Map(amap.data, amap.meta)
        assert isinstance(pair_map, sunpy.map.GenericMap)
         
    
    def test_save(self):
        #Test save out
        eitmap = sunpy.Map(a_fname)
        eitmap.save("eit_save.fits", filetype='fits', clobber=True)
        backin = sunpy.Map("eit_save.fits")
        assert isinstance(backin, sunpy.map.sources.EITMap)
        os.remove("eit_save.fits")
    
#==============================================================================
# Sources Tests
#==============================================================================
    def test_sdo(self):
        #Test an AIAMap
        aia = sunpy.Map(sunpy.AIA_171_IMAGE)
        assert isinstance(aia,sunpy.map.sources.AIAMap)
        #Test a HMIMap
        
    def test_soho(self):
        #Test EITMap, LASCOMap & MDIMap
        eit = sunpy.Map(filepath + "/EIT/efz20040301.000010_s.fits")
        assert isinstance(eit,sunpy.map.sources.EITMap)
    
        lasco = sunpy.Map(filepath + "/lasco_c2_25299383_s.fts")
        assert isinstance(lasco,sunpy.map.sources.LASCOMap)
        
        mdi_c = sunpy.Map(filepath + "/mdi_fd_Ic_6h_01d.5871.0000_s.fits")
        assert isinstance(mdi_c,sunpy.map.sources.MDIMap)
        
        mdi_m = sunpy.Map(filepath + "/mdi_fd_M_96m_01d.5874.0005_s.fits")
        assert isinstance(mdi_m,sunpy.map.sources.MDIMap)
        
    def test_stereo(self):    
        #Test EUVIMap & CORMap
        euvi = sunpy.Map(filepath + "/euvi_20090615_000900_n4euA_s.fts")
        assert isinstance(euvi,sunpy.map.sources.EUVIMap)
        
        cor = sunpy.Map(filepath + "/cor1_20090615_000500_s4c1A.fts")
        assert isinstance(cor,sunpy.map.sources.CORMap)
        
    def test_rhessi(self):    
        #Test RHESSIMap
        rhessi = sunpy.Map(sunpy.RHESSI_IMAGE)
        assert isinstance(rhessi,sunpy.map.sources.RHESSIMap)
        
        #Test SWAPMap
        
        #Test XRTMap
        
        #Test SXTMap