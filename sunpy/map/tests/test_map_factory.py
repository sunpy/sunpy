# -*- coding: utf-8 -*-
"""
Created on Fri Jun 21 15:05:09 2013

@author: stuart
"""
import os
import glob
import numpy as np

import pytest

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
        cube = sunpy.map.Map(a_list_of_many, cube=True)
        assert isinstance(cube, sunpy.map.MapCube)

    def test_composite(self):
        #Test making a CompositeMap
        comp = sunpy.map.Map(sunpy.AIA_171_IMAGE, sunpy.RHESSI_IMAGE,
                         composite=True)
        assert isinstance(comp, sunpy.map.CompositeMap)

    def test_patterns(self):
        ## Test different Map pattern matching ##
        # File name
        eitmap = sunpy.map.Map(a_fname)
        assert isinstance(eitmap, sunpy.map.GenericMap)
        # Directory
        maps = sunpy.map.Map(os.path.join(filepath,"EIT"))
        assert isinstance(maps, list)
        assert ([isinstance(amap,sunpy.map.GenericMap) for amap in maps])
        # Glob
        maps = sunpy.map.Map(os.path.join(filepath,"EIT")+'/*')
        assert isinstance(maps, list)
        assert ([isinstance(amap,sunpy.map.GenericMap) for amap in maps])
        # Already a Map
        amap = sunpy.map.Map(maps[0])
        assert isinstance(amap, sunpy.map.GenericMap)
        # A list of filenames
        maps = sunpy.map.Map(a_list_of_many)
        assert isinstance(maps, list)
        assert ([isinstance(amap,sunpy.map.GenericMap) for amap in maps])
        # Data-header pair in a tuple
        pair_map = sunpy.map.Map((amap.data, amap.meta))
        assert isinstance(pair_map, sunpy.map.GenericMap)
        # Data-header pair not in a tuple
        pair_map = sunpy.map.Map(amap.data, amap.meta)
        assert isinstance(pair_map, sunpy.map.GenericMap)
        #Custom Map
        data = np.arange(0,100).reshape(10,10)
        header = {'cdelt1': 10, 'cdelt2': 10, 'telescop':'sunpy'}
        pair_map = sunpy.map.Map(data, header)
        assert isinstance(pair_map, sunpy.map.GenericMap)

    def test_databaseentry(self):
        db = sunpy.database.Database(url='sqlite://', default_waveunit='angstrom')
        db.add_from_file(a_fname)

        res = db.get_entry_by_id(1)
        db_map = sunpy.map.Map(res)
        assert isinstance(db_map, sunpy.map.GenericMap)

    @pytest.mark.online
    def test_url_pattern(self):
        # A URL
        amap = sunpy.map.Map("https://raw.github.com/sunpy/sunpy/master/sunpy/data/sample/AIA20110319_105400_0171.fits")
        assert isinstance(amap, sunpy.map.GenericMap)

    def test_save(self):
        #Test save out
        eitmap = sunpy.map.Map(a_fname)
        eitmap.save("eit_save.fits", filetype='fits', clobber=True)
        backin = sunpy.map.Map("eit_save.fits")
        assert isinstance(backin, sunpy.map.sources.EITMap)
        os.remove("eit_save.fits")

#==============================================================================
# Sources Tests
#==============================================================================
    def test_sdo(self):
        #Test an AIAMap
        aia = sunpy.map.Map(sunpy.AIA_171_IMAGE)
        assert isinstance(aia,sunpy.map.sources.AIAMap)
        #Test a HMIMap

    def test_soho(self):
        #Test EITMap, LASCOMap & MDIMap
        eit = sunpy.map.Map(filepath + "/EIT/efz20040301.000010_s.fits")
        assert isinstance(eit,sunpy.map.sources.EITMap)

        lasco = sunpy.map.Map(filepath + "/lasco_c2_25299383_s.fts")
        assert isinstance(lasco,sunpy.map.sources.LASCOMap)

        mdi_c = sunpy.map.Map(filepath + "/mdi_fd_Ic_6h_01d.5871.0000_s.fits")
        assert isinstance(mdi_c,sunpy.map.sources.MDIMap)

        mdi_m = sunpy.map.Map(filepath + "/mdi_fd_M_96m_01d.5874.0005_s.fits")
        assert isinstance(mdi_m,sunpy.map.sources.MDIMap)

    def test_stereo(self):
        #Test EUVIMap & CORMap
        euvi = sunpy.map.Map(filepath + "/euvi_20090615_000900_n4euA_s.fts")
        assert isinstance(euvi,sunpy.map.sources.EUVIMap)

        cor = sunpy.map.Map(filepath + "/cor1_20090615_000500_s4c1A.fts")
        assert isinstance(cor,sunpy.map.sources.CORMap)

    def test_rhessi(self):
        #Test RHESSIMap
        rhessi = sunpy.map.Map(sunpy.RHESSI_IMAGE)
        assert isinstance(rhessi,sunpy.map.sources.RHESSIMap)

        #Test SWAPMap

        #Test XRTMap

        #Test SXTMap
