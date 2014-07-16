# -*- coding: utf-8 -*-
import numpy as np
import os

import pytest

import sunpy
import sunpy.io
import sunpy.data.test

# SunPy's JPEG2000 capabilities rely on the glymur library.  First we check to
# make sure that glymur imports correctly before proceeding.
try:
    import glymur
except ImportError:
    glymur_imports = False
else:
    glymur_imports = True

# Skip ana tests if we are on Windows or we can't import the c extension.
import platform
if platform.system() == 'Windows':
    skip_ana = True
else:
    skip_ana = False

try:
    import sunpy.io._pyana
except ImportError:
    skip_ana = True
else:
    skip_ana = skip_ana or False
    
#==============================================================================
# Test, read, get_header and write through the file independant layer
#==============================================================================
class TestFiletools():

    def test_read_file_fits(self):
        #Test read FITS
        aiapair = sunpy.io.read_file(sunpy.AIA_171_IMAGE)
        assert isinstance(aiapair, list)
        assert len(aiapair) == 1
        assert len(aiapair[0]) == 2
        assert isinstance(aiapair[0][0], np.ndarray)
        assert isinstance(aiapair[0][1], sunpy.io.header.FileHeader)

        #Test read multi HDU list
        pairs = sunpy.io.read_file(sunpy.RHESSI_IMAGE)
        assert isinstance(pairs, list)
        assert len(pairs) == 4
        assert all([len(p) == 2 for p in pairs])
        assert all([isinstance(p[0], np.ndarray) for p in pairs])
        assert all([isinstance(p[1],
                               sunpy.io.header.FileHeader) for p in pairs])

    @pytest.mark.skipif("glymur_imports is False")
    def test_read_file_jp2(self):
        #Test read jp2
        pair = sunpy.io.read_file(os.path.join(sunpy.data.test.rootdir,
                               "2013_06_24__17_31_30_84__SDO_AIA_AIA_193.jp2"))

        assert isinstance(pair, list)
        assert len(pair) == 1
        assert len(pair[0]) == 2
        assert isinstance(pair[0][0], np.ndarray)
        assert isinstance(pair[0][1], sunpy.io.header.FileHeader)

    def test_read_file_header_fits(self):
        #Test FITS
        hlist = sunpy.io.read_file_header(sunpy.AIA_171_IMAGE)
        assert isinstance(hlist, list)
        assert len(hlist) == 1
        assert isinstance(hlist[0], sunpy.io.header.FileHeader)

    @pytest.mark.skipif("glymur_imports is False")
    def test_read_file_header_jp2(self):
        #Test jp2
        hlist = sunpy.io.read_file_header(os.path.join(sunpy.data.test.rootdir,
                               "2013_06_24__17_31_30_84__SDO_AIA_AIA_193.jp2"))
        assert isinstance(hlist, list)
        assert len(hlist) == 1
        assert isinstance(hlist[0], sunpy.io.header.FileHeader)


    def test_write_file_fits(self):
        #Test write FITS
        aiapair = sunpy.io.read_file(sunpy.AIA_171_IMAGE)[0]
        sunpy.io.write_file("aia_171_image.fits", aiapair[0], aiapair[1],
                            clobber=True)
        assert os.path.exists("aia_171_image.fits")
        outpair = sunpy.io.read_file(sunpy.AIA_171_IMAGE)[0]
        assert np.all(np.equal(outpair[0], aiapair[0]))
        assert outpair[1] == aiapair[1]
        os.remove("aia_171_image.fits")

    @pytest.mark.skipif("skip_ana is True")
    def test_read_file_ana(self):
        ana_data = sunpy.io.read_file(os.path.join(sunpy.data.test.rootdir,"test_ana.fz"))
        assert isinstance(ana_data, list)
        assert len(ana_data) == 1
        assert len(ana_data[0]) == 2
        assert isinstance(ana_data[0][0], np.ndarray)
        assert isinstance(ana_data[0][1], sunpy.io.header.FileHeader)

    @pytest.mark.skipif("skip_ana is True")
    def test_read_file__header_ana(self):
        ana_data = sunpy.io.read_file_header(os.path.join(sunpy.data.test.rootdir,"test_ana.fz"))
        assert isinstance(ana_data, list)
        assert len(ana_data) == 1
        assert isinstance(ana_data[0], sunpy.io.header.FileHeader)

    @pytest.mark.skipif("skip_ana is True")
    def test_write_file_ana(self):
        ana = sunpy.io.read_file(os.path.join(sunpy.data.test.rootdir,"test_ana.fz"))[0]
        sunpy.io.write_file("ana_test_write.fz", ana[0], str(ana[1]))
        assert os.path.exists("ana_test_write.fz")
        outpair = sunpy.io.read_file(os.path.join(sunpy.data.test.rootdir,"test_ana.fz"))
        assert np.all(np.equal(outpair[0][1], ana[1]))
        assert outpair[0][1] == ana[1]
        os.remove("ana_test_write.fz")

    #TODO: Test write jp2
