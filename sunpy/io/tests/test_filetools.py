# -*- coding: utf-8 -*-
import numpy as np
import os
<<<<<<< HEAD
from pathlib import Path
=======
import pathlib
>>>>>>> pathlib added astropy files not touched yet

import sunpy
import sunpy.io
import sunpy.data.test

from sunpy.tests.helpers import skip_glymur, skip_ana

testpath = sunpy.data.test.rootdir

<<<<<<< HEAD
RHESSI_IMAGE = str(Path.home().joinpath(testpath, 'hsi_image_20101016_191218.fits'))
EIT_195_IMAGE = str(Path.home().joinpath(testpath, 'EIT/efz20040301.000010_s.fits'))
AIA_171_IMAGE = str(Path.home().joinpath(testpath, 'aia_171_level1.fits'))
=======
RHESSI_IMAGE = str(pathlib.Path.home().joinpath(testpath, 'hsi_image_20101016_191218.fits'))
EIT_195_IMAGE = str(pathlib.Path.home().joinpath(testpath, 'EIT/efz20040301.000010_s.fits'))
AIA_171_IMAGE = str(pathlib.Path.home().joinpath(testpath, 'aia_171_level1.fits'))
>>>>>>> pathlib added astropy files not touched yet

#==============================================================================
# Test, read, get_header and write through the file independent layer
#==============================================================================
class TestFiletools(object):

    def test_read_file_fits(self):
        #Test read FITS
        aiapair = sunpy.io.read_file(AIA_171_IMAGE)
        assert isinstance(aiapair, list)
        assert len(aiapair) == 1
        assert len(aiapair[0]) == 2
        assert isinstance(aiapair[0][0], np.ndarray)
        assert isinstance(aiapair[0][1], sunpy.io.header.FileHeader)

        #Test read multi HDU list
        pairs = sunpy.io.read_file(RHESSI_IMAGE)
        assert isinstance(pairs, list)
        assert len(pairs) == 4
        assert all([len(p) == 2 for p in pairs])
        assert all([isinstance(p[0], np.ndarray) for p in pairs])
        assert all([isinstance(p[1],
                               sunpy.io.header.FileHeader) for p in pairs])

    def test_read_file_fits_gzip(self):
        # Test read gzipped fits file
        for fits_extension in [".fts", ".fit", ".fits"]:
<<<<<<< HEAD
            pair = sunpy.io.read_file(str(Path.home().joinpath(sunpy.data.test.rootdir, "gzip_test{ext}.gz".format(ext=fits_extension))))
=======
            pair = sunpy.io.read_file(str(pathlib.Path.home().joinpath(sunpy.data.test.rootdir, "gzip_test{ext}.gz".format(ext=fits_extension))))
>>>>>>> pathlib added astropy files not touched yet
            assert isinstance(pair, list)
            assert len(pair) == 1
            assert len(pair[0]) == 2
            assert isinstance(pair[0][0], np.ndarray)
            assert isinstance(pair[0][1], sunpy.io.header.FileHeader)
            assert np.all(pair[0][0] == np.tile(np.arange(32), (32, 1)).transpose())

    @skip_glymur
    def test_read_file_jp2(self):
        #Test read jp2
<<<<<<< HEAD
        pair = sunpy.io.read_file(str(Path.home().joinpath(sunpy.data.test.rootdir,
=======
        pair = sunpy.io.read_file(str(pathlib.Path.home().joinpath(sunpy.data.test.rootdir,
>>>>>>> pathlib added astropy files not touched yet
                               "2013_06_24__17_31_30_84__SDO_AIA_AIA_193.jp2")))

        assert isinstance(pair, list)
        assert len(pair) == 1
        assert len(pair[0]) == 2
        assert isinstance(pair[0][0], np.ndarray)
        assert isinstance(pair[0][1], sunpy.io.header.FileHeader)

    def test_read_file_header_fits(self):
        # Test FITS
        hlist = sunpy.io.read_file_header(AIA_171_IMAGE)
        assert isinstance(hlist, list)
        assert len(hlist) == 1
        assert isinstance(hlist[0], sunpy.io.header.FileHeader)

    @skip_glymur
    def test_read_file_header_jp2(self):
        # Test jp2
<<<<<<< HEAD
        hlist = sunpy.io.read_file_header(str(Path.home().joinpath(sunpy.data.test.rootdir,
=======
        hlist = sunpy.io.read_file_header(str(pathlib.Path.home().joinpath(sunpy.data.test.rootdir,
>>>>>>> pathlib added astropy files not touched yet
                                            "2013_06_24__17_31_30_84__SDO_AIA_AIA_193.jp2")))
        assert isinstance(hlist, list)
        assert len(hlist) == 1
        assert isinstance(hlist[0], sunpy.io.header.FileHeader)

    def test_write_file_fits(self):
        # Test write FITS
        aiapair = sunpy.io.read_file(AIA_171_IMAGE)[0]
        sunpy.io.write_file("aia_171_image.fits", aiapair[0], aiapair[1],
                            clobber=True)
<<<<<<< HEAD
        assert Path("aia_171_image.fits").exists()
=======
        assert pathlib.Path("aia_171_image.fits").exists
>>>>>>> pathlib added astropy files not touched yet
        outpair = sunpy.io.read_file(AIA_171_IMAGE)[0]
        assert np.all(np.equal(outpair[0], aiapair[0]))
        assert outpair[1] == aiapair[1]
        os.remove("aia_171_image.fits")

    def test_write_file_fits_bytes(self):
        # Test write FITS
        aiapair = sunpy.io.read_file(AIA_171_IMAGE)[0]
        with open("aia_171_image_bytes.fits", 'wb') as fileo:
            sunpy.io.write_file(fileo, aiapair[0], aiapair[1], filetype='fits')
<<<<<<< HEAD
        assert Path("aia_171_image_bytes.fits").exists
=======
        assert pathlib.Path("aia_171_image_bytes.fits").exists
>>>>>>> pathlib added astropy files not touched yet
        outpair = sunpy.io.read_file(AIA_171_IMAGE)[0]
        assert np.all(np.equal(outpair[0], aiapair[0]))
        assert outpair[1] == aiapair[1]
        os.remove("aia_171_image_bytes.fits")

    @skip_ana
    def test_read_file_ana(self):
<<<<<<< HEAD
        ana_data = sunpy.io.read_file(str(Path.home().joinpath(sunpy.data.test.rootdir,"test_ana.fz")))
=======
        ana_data = sunpy.io.read_file(str(pathlib.Path.home().joinpath(sunpy.data.test.rootdir,"test_ana.fz")))
>>>>>>> pathlib added astropy files not touched yet
        assert isinstance(ana_data, list)
        assert len(ana_data) == 1
        assert len(ana_data[0]) == 2
        assert isinstance(ana_data[0][0], np.ndarray)
        assert isinstance(ana_data[0][1], sunpy.io.header.FileHeader)

    @skip_ana
    def test_read_file__header_ana(self):
<<<<<<< HEAD
        ana_data = sunpy.io.read_file_header(str(Path.home().joinpath(sunpy.data.test.rootdir,"test_ana.fz")))
=======
        ana_data = sunpy.io.read_file_header(str(pathlib.Path.home().joinpath(sunpy.data.test.rootdir,"test_ana.fz")))
>>>>>>> pathlib added astropy files not touched yet
        assert isinstance(ana_data, list)
        assert len(ana_data) == 1
        assert isinstance(ana_data[0], sunpy.io.header.FileHeader)

    @skip_ana
    def test_write_file_ana(self):
<<<<<<< HEAD
        ana = sunpy.io.read_file(str(Path.home().joinpath(sunpy.data.test.rootdir,"test_ana.fz")))[0]
        sunpy.io.write_file("ana_test_write.fz", ana[0], str(ana[1]))
        assert Path("ana_test_write.fz").exists
        outpair = sunpy.io.read_file(str(Path.home().joinpath(sunpy.data.test.rootdir,"test_ana.fz")))
=======
        ana = sunpy.io.read_file(str(pathlib.Path.home().joinpath(sunpy.data.test.rootdir,"test_ana.fz")))[0]
        sunpy.io.write_file("ana_test_write.fz", ana[0], str(ana[1]))
        assert pathlib.Path("ana_test_write.fz").exists
        outpair = sunpy.io.read_file(str(pathlib.Path.home().joinpath(sunpy.data.test.rootdir,"test_ana.fz")))
>>>>>>> pathlib added astropy files not touched yet
        assert np.all(np.equal(outpair[0][1], ana[1]))
        assert outpair[0][1] == ana[1]
        os.remove("ana_test_write.fz")

    #TODO: Test write jp2
