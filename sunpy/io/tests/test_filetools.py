import os
import pathlib

import numpy as np
import pytest

import sunpy
import sunpy.io
from sunpy.data.test import get_test_filepath
from sunpy.tests.helpers import skip_ana, skip_glymur

TEST_RHESSI_IMAGE = get_test_filepath('hsi_image_20101016_191218.fits')
TEST_AIA_IMAGE = get_test_filepath('aia_171_level1.fits')
# Some of the tests images contain an invalid BLANK keyword;
pytestmark = pytest.mark.filterwarnings("ignore:Invalid 'BLANK' keyword in header")


def test_read_file_network_fits(mocker):
    # Aim is to verify that we can read a files from a URL
    # But it is mocked to prevent network access
    url = "https://hesperia.gsfc.nasa.gov/rhessi_extras/imagecube_fits/2015/12/20/20151220_2228_2248/hsi_imagecube_clean_20151220_2228_13tx3e.fits"
    mock = mocker.patch("astropy.io.fits.file.download_file", return_value=TEST_AIA_IMAGE)
    data = sunpy.io.read_file(url)
    assert mock.call_args[0] == (url,)
    assert isinstance(data, list)
    assert len(data) == 1
    assert len(data[0]) == 2
    assert isinstance(data[0][0], np.ndarray)
    assert isinstance(data[0][1], sunpy.io.header.FileHeader)


def test_read_file_fits():
    # Aim is to verify that we can read a FITS file
    aiapair = sunpy.io.read_file(TEST_AIA_IMAGE)
    assert isinstance(aiapair, list)
    assert len(aiapair) == 1
    assert len(aiapair[0]) == 2
    assert isinstance(aiapair[0][0], np.ndarray)
    assert isinstance(aiapair[0][1], sunpy.io.header.FileHeader)


def test_read_file_fits_multple_hdu():
    # Aim is to verify that we can read a FITS file with multiple HDUs
    pairs = sunpy.io.read_file(TEST_RHESSI_IMAGE)
    assert isinstance(pairs, list)
    assert len(pairs) == 4
    assert all(len(p) == 2 for p in pairs)
    assert all(isinstance(p[0], np.ndarray) for p in pairs)
    assert all(isinstance(p[1], sunpy.io.header.FileHeader) for p in pairs)


@pytest.mark.parametrize('fname', ["gzip_test.fts.gz", "gzip_test.fits.gz", "gzip_test.fit.gz", "gzip_fits_test.file"])
def test_read_file_fits_gzip(fname):
    # Aim is to verify that we from a gzipped FITS file
    # that has a range of different file extensions
    pair = sunpy.io.read_file(get_test_filepath(fname))
    assert isinstance(pair, list)
    assert len(pair) == 1
    assert len(pair[0]) == 2
    assert isinstance(pair[0][0], np.ndarray)
    assert isinstance(pair[0][1], sunpy.io.header.FileHeader)
    assert np.all(pair[0][0] == np.tile(np.arange(32), (32, 1)).transpose())


@skip_glymur
def test_read_file_jp2():
    # Aim is to verify that we can read a JP2 file
    pair = sunpy.io.read_file(get_test_filepath("2013_06_24__17_31_30_84__SDO_AIA_AIA_193.jp2"))
    assert isinstance(pair, list)
    assert len(pair) == 1
    assert len(pair[0]) == 2
    assert isinstance(pair[0][0], np.ndarray)
    assert isinstance(pair[0][1], sunpy.io.header.FileHeader)


@skip_glymur
def test_read_file_header_jp2():
    # Aim is to verify that we can read a header from a JP2 file
    hlist = sunpy.io.read_file_header(get_test_filepath("2013_06_24__17_31_30_84__SDO_AIA_AIA_193.jp2"))
    assert isinstance(hlist, list)
    assert len(hlist) == 1
    assert isinstance(hlist[0], sunpy.io.header.FileHeader)


def test_read_file_header_fits():
    # Aim is to verify that we can read a header from a FITS file
    hlist = sunpy.io.read_file_header(TEST_AIA_IMAGE)
    assert isinstance(hlist, list)
    assert len(hlist) == 1
    assert isinstance(hlist[0], sunpy.io.header.FileHeader)


@skip_ana
def test_read_file_ana():
    # Aim is to verify that we can read a ANA file
    ana_data = sunpy.io.read_file(get_test_filepath("test_ana.fz"))
    assert isinstance(ana_data, list)
    assert len(ana_data) == 1
    assert len(ana_data[0]) == 2
    assert isinstance(ana_data[0][0], np.ndarray)
    assert isinstance(ana_data[0][1], sunpy.io.header.FileHeader)


@skip_ana
def test_read_file_header_ana():
    # Aim is to verify that we can read a header from a ANA file
    ana_data = sunpy.io.read_file_header(get_test_filepath("test_ana.fz"))
    assert isinstance(ana_data, list)
    assert len(ana_data) == 1
    assert isinstance(ana_data[0], sunpy.io.header.FileHeader)


@skip_ana
def test_write_file_ana():
    # Aim is to verify that we can write a ANA file and read back correctly
    ana = sunpy.io.read_file(get_test_filepath("test_ana.fz"))[0]
    sunpy.io.write_file("ana_test_write.fz", ana[0], str(ana[1]))
    assert os.path.exists("ana_test_write.fz")
    outpair = sunpy.io.read_file(get_test_filepath("test_ana.fz"))
    assert np.all(np.equal(outpair[0][1], ana[1]))
    assert outpair[0][1] == ana[1]


@pytest.mark.parametrize('fname', ['aia_171_image.fits', pathlib.Path('aia_171_image.fits')])
def test_write_file_fits(fname, tmpdir):
    # Aim is to verify that we can write a FITS file and read it back correctly
    aiapair = sunpy.io.read_file(TEST_AIA_IMAGE)[0]
    filepath = tmpdir / fname
    sunpy.io.write_file(filepath, aiapair[0], aiapair[1])
    assert filepath.exists()
    outpair = sunpy.io.read_file(filepath)[0]
    assert np.all(np.equal(outpair[0], aiapair[0]))
    assert outpair[1] == aiapair[1]


def test_write_file_fits_bytes(tmpdir):
    # Aim is to verify that we can write a FITS file via a file-like object and read it back correctly
    aiapair = sunpy.io.read_file(TEST_AIA_IMAGE)[0]
    filepath = tmpdir / "aia_171_image_bytes.fits"
    with open(filepath, "wb") as fileo:
        sunpy.io.write_file(fileo, aiapair[0], aiapair[1], filetype='fits')
    assert filepath.exists()
    outpair = sunpy.io.read_file(filepath)[0]
    assert np.all(np.equal(outpair[0], aiapair[0]))
    assert outpair[1] == aiapair[1]

def test_missing_file_extension(tmpdir):
    # Aim is read a FITS file without an extension
    aiapair = sunpy.io.read_file(TEST_AIA_IMAGE)[0]
    filepath = tmpdir / "test"
    with open(filepath, 'wb') as fileo:
        sunpy.io.write_file(fileo, aiapair[0], aiapair[1], filetype='fits')
    assert filepath.exists()
    outpair = sunpy.io.read_file(filepath)[0]
    assert np.all(np.equal(outpair[0], aiapair[0]))
    assert outpair[1] == aiapair[1]
