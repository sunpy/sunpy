import os
import pathlib

import numpy as np
import pytest

import sunpy
import sunpy.data.test
import sunpy.io
from sunpy.tests.helpers import skip_ana, skip_glymur

testpath = sunpy.data.test.rootdir

RHESSI_IMAGE = os.path.join(testpath, 'hsi_image_20101016_191218.fits')
EIT_195_IMAGE = os.path.join(testpath, 'EIT/efz20040301.000010_s.fits')
AIA_171_IMAGE = os.path.join(testpath, 'aia_171_level1.fits')


# Some of the tests iamges contain an invalid BLANK keyword; ignore the warning
# raised by this
pytestmark = pytest.mark.filterwarnings("ignore:Invalid 'BLANK' keyword in header")


# =============================================================================
# Test, read, get_header and write through the file independent layer
# =============================================================================
class TestFiletools:

    def test_read_file_fits(self):
        # Test read FITS
        aiapair = sunpy.io.read_file(AIA_171_IMAGE)
        assert isinstance(aiapair, list)
        assert len(aiapair) == 1
        assert len(aiapair[0]) == 2
        assert isinstance(aiapair[0][0], np.ndarray)
        assert isinstance(aiapair[0][1], sunpy.io.header.FileHeader)

        # Test read multi HDU list
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
            pair = sunpy.io.read_file(os.path.join(
                sunpy.data.test.rootdir, f"gzip_test{fits_extension}.gz"))
            assert isinstance(pair, list)
            assert len(pair) == 1
            assert len(pair[0]) == 2
            assert isinstance(pair[0][0], np.ndarray)
            assert isinstance(pair[0][1], sunpy.io.header.FileHeader)
            assert np.all(pair[0][0] == np.tile(np.arange(32), (32, 1)).transpose())

    @skip_glymur
    def test_read_file_jp2(self):
        # Test read jp2
        pair = sunpy.io.read_file(os.path.join(sunpy.data.test.rootdir,
                                               "2013_06_24__17_31_30_84__SDO_AIA_AIA_193.jp2"))

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
        hlist = sunpy.io.read_file_header(os.path.join(sunpy.data.test.rootdir,
                                                       "2013_06_24__17_31_30_84__SDO_AIA_AIA_193.jp2"))
        assert isinstance(hlist, list)
        assert len(hlist) == 1
        assert isinstance(hlist[0], sunpy.io.header.FileHeader)

    @pytest.mark.parametrize('fname', ['aia_171_image.fits',
                                       pathlib.Path('aia_171_image.fits')])
    def test_write_file_fits(self, fname):
        # Test write FITS
        aiapair = sunpy.io.read_file(AIA_171_IMAGE)[0]
        sunpy.io.write_file(fname, aiapair[0], aiapair[1],
                            overwrite=True)
        assert os.path.exists("aia_171_image.fits")
        outpair = sunpy.io.read_file(AIA_171_IMAGE)[0]
        assert np.all(np.equal(outpair[0], aiapair[0]))
        assert outpair[1] == aiapair[1]
        os.remove("aia_171_image.fits")

    def test_write_file_fits_bytes(self):
        # Test write FITS
        aiapair = sunpy.io.read_file(AIA_171_IMAGE)[0]
        with open("aia_171_image_bytes.fits", 'wb') as fileo:
            sunpy.io.write_file(fileo, aiapair[0], aiapair[1], filetype='fits')
        assert os.path.exists("aia_171_image_bytes.fits")
        outpair = sunpy.io.read_file(AIA_171_IMAGE)[0]
        assert np.all(np.equal(outpair[0], aiapair[0]))
        assert outpair[1] == aiapair[1]
        os.remove("aia_171_image_bytes.fits")

    @skip_ana
    def test_read_file_ana(self):
        ana_data = sunpy.io.read_file(os.path.join(sunpy.data.test.rootdir, "test_ana.fz"))
        assert isinstance(ana_data, list)
        assert len(ana_data) == 1
        assert len(ana_data[0]) == 2
        assert isinstance(ana_data[0][0], np.ndarray)
        assert isinstance(ana_data[0][1], sunpy.io.header.FileHeader)

    @skip_ana
    def test_read_file__header_ana(self):
        ana_data = sunpy.io.read_file_header(os.path.join(sunpy.data.test.rootdir, "test_ana.fz"))
        assert isinstance(ana_data, list)
        assert len(ana_data) == 1
        assert isinstance(ana_data[0], sunpy.io.header.FileHeader)

    @skip_ana
    def test_write_file_ana(self):
        ana = sunpy.io.read_file(os.path.join(sunpy.data.test.rootdir, "test_ana.fz"))[0]
        sunpy.io.write_file("ana_test_write.fz", ana[0], str(ana[1]))
        assert os.path.exists("ana_test_write.fz")
        outpair = sunpy.io.read_file(os.path.join(sunpy.data.test.rootdir, "test_ana.fz"))
        assert np.all(np.equal(outpair[0][1], ana[1]))
        assert outpair[0][1] == ana[1]
        os.remove("ana_test_write.fz")

    # TODO: Test write jp2
