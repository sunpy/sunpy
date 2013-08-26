# -*- coding: utf-8 -*-
import numpy as np
import os

import sunpy
import sunpy.io
import sunpy.data.test
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

    def test_read_file_jp2(self):
        #Test read jp2
        pair = sunpy.io.read_file(os.path.join(sunpy.data.test.rootdir,
                               "2013_06_24__17_31_30_84__SDO_AIA_AIA_193.jp2"),
                               j2k_to_image='j2k_to_image')
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

    #Test write jp2
