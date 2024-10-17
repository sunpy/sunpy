import os
import pathlib
from unittest.mock import patch

import numpy as np
import pytest

from sunpy.data.test import get_test_filepath
from sunpy.io._file_tools import read_file, read_file_header, write_file
from sunpy.io._header import FileHeader
from sunpy.tests.helpers import skip_ana, skip_glymur

TEST_RHESSI_IMAGE = get_test_filepath('hsi_image_20101016_191218.fits')
TEST_AIA_IMAGE = get_test_filepath('aia_171_level1.fits')
TEST_AIA_IMAGE_URI = "s3://data.sunpy.org/sunpy/AIA20110607_065843_0193_cutout.fits"

# Some of the tests images contain an invalid BLANK keyword;
pytestmark = pytest.mark.filterwarnings("ignore:Invalid 'BLANK' keyword in header")


def test_read_file_network_fits():
    # Aim is to verify that we can read a files from a URL
    # but it is mocked to prevent network access
    url = "https://hesperia.gsfc.nasa.gov/rhessi_extras/imagecube_fits/2015/12/20/20151220_2228_2248/hsi_imagecube_clean_20151220_2228_13tx3e.fits"
    with patch("astropy.io.fits.file.download_file", return_value=TEST_AIA_IMAGE) as mock:
        hdulist = read_file(url)
        assert mock.call_args[0][0] == url
    assert isinstance(hdulist, list)
    assert len(hdulist) == 1
    assert len(hdulist[0]) == 2
    assert isinstance(hdulist[0][0], np.ndarray)
    assert isinstance(hdulist[0][1], FileHeader)


def test_read_file_fits():
    # Aim is to verify that we can read a FITS file
    aia_header, aia_data = read_file(TEST_AIA_IMAGE)[0]
    assert isinstance(aia_header, np.ndarray)
    assert isinstance(aia_data, FileHeader)

@pytest.mark.remote_data
def test_read_file_fits_uri():
    # Aim is to verify that we can read a FITS file from s3
    aia_header, aia_data = read_file(TEST_AIA_IMAGE_URI, fsspec_kwargs={'anon':True})[0]
    assert isinstance(aia_header, np.ndarray)
    assert isinstance(aia_data, FileHeader)

def test_read_file_fits_multiple_hdu():
    # Aim is to verify that we can read a FITS file with multiple HDUs
    hdulist = read_file(TEST_RHESSI_IMAGE)
    assert isinstance(hdulist, list)
    assert len(hdulist) == 4
    assert all(len(hdupair) == 2 for hdupair in hdulist)
    assert all(isinstance(hdupair[0], np.ndarray) for hdupair in hdulist)
    assert all(isinstance(hdupair[1], FileHeader) for hdupair in hdulist)


@pytest.mark.parametrize('fname',
                         ["gzip_test.fts.gz", "gzip_test.fits.gz", "gzip_test.fit.gz", "gzip_fits_test.file"])
def test_read_file_fits_gzip(fname):
    # Aim is to verify that we from a gzipped FITS file
    # that has a range of different file extensions
    hdulist = read_file(get_test_filepath(fname))
    assert isinstance(hdulist, list)
    assert len(hdulist) == 1
    assert len(hdulist[0]) == 2
    assert isinstance(hdulist[0][0], np.ndarray)
    assert isinstance(hdulist[0][1], FileHeader)
    assert np.all(hdulist[0][0] == np.tile(np.arange(32), (32, 1)).transpose())


def test_read_file_fits_fz_ext(tmp_path):
    test_fits = tmp_path / 'test.fits.fz'
    with open(TEST_AIA_IMAGE, 'rb') as orig:
        with (test_fits).open('wb') as test:
            test.write(orig.read())

    aia_header, aia_data = read_file(str(test_fits))[0]
    assert isinstance(aia_header, np.ndarray)
    assert isinstance(aia_data, FileHeader)

@skip_glymur
def test_read_file_jp2():
    # Aim is to verify that we can read a JP2 file
    hdulist = read_file(get_test_filepath("2013_06_24__17_31_30_84__SDO_AIA_AIA_193.jp2"))
    assert isinstance(hdulist, list)
    assert len(hdulist) == 1
    assert len(hdulist[0]) == 2
    assert isinstance(hdulist[0][0], np.ndarray)
    assert isinstance(hdulist[0][1], FileHeader)


@skip_glymur
def test_read_file_header_jp2():
    # Aim is to verify that we can read a header from a JP2 file
    hdulist = read_file_header(get_test_filepath("2013_06_24__17_31_30_84__SDO_AIA_AIA_193.jp2"))
    assert isinstance(hdulist, list)
    assert len(hdulist) == 1
    assert isinstance(hdulist[0], FileHeader)


def test_read_file_header_fits():
    # Aim is to verify that we can read a header from a FITS file
    hdulist = read_file_header(TEST_AIA_IMAGE)
    assert isinstance(hdulist, list)
    assert len(hdulist) == 1
    assert isinstance(hdulist[0], FileHeader)


@skip_ana
@pytest.mark.filterwarnings("ignore::sunpy.util.exceptions.SunpyDeprecationWarning")
def test_read_file_ana():
    # Aim is to verify that we can read a ANA file
    hdulist = read_file(get_test_filepath("test_ana.fz"))
    assert isinstance(hdulist, list)
    assert len(hdulist) == 1
    assert len(hdulist[0]) == 2
    assert isinstance(hdulist[0][0], np.ndarray)
    assert isinstance(hdulist[0][1], FileHeader)


@skip_ana
@pytest.mark.filterwarnings("ignore::sunpy.util.exceptions.SunpyDeprecationWarning")
def test_read_file_header_ana():
    # Aim is to verify that we can read a header from a ANA file
    hdulist = read_file_header(get_test_filepath("test_ana.fz"))
    assert isinstance(hdulist, list)
    assert len(hdulist) == 1
    assert isinstance(hdulist[0], FileHeader)


@skip_ana
@pytest.mark.filterwarnings("ignore::sunpy.util.exceptions.SunpyDeprecationWarning")
def test_write_file_ana(tmpdir):
    # Aim is to verify that we can write a ANA file and read back correctly
    ana_header, ana_data = read_file(get_test_filepath("test_ana.fz"))[0][::-1]
    write_file(str(tmpdir.join("ana_test_write.fz")), ana_data, str(ana_header))
    assert os.path.exists(str(tmpdir.join("ana_test_write.fz")))
    test_ana_header, test_ana_data = read_file(get_test_filepath("test_ana.fz"))[0][::-1]
    assert np.all(np.equal(test_ana_data, ana_data))
    assert test_ana_header == ana_header


@pytest.mark.parametrize('fname',
                         ['aia_171_image.fits', pathlib.Path('aia_171_image.fits')])
def test_write_file_fits(fname, tmpdir):
    # Aim is to verify that we can write a FITS file and read it back correctly
    aia_header, aia_data = read_file(TEST_AIA_IMAGE)[0][::-1]
    filepath = tmpdir / fname
    write_file(filepath, aia_data, aia_header)
    assert filepath.exists()
    test_aia_header, test_aia_data = read_file(filepath)[0][::-1]
    assert np.all(np.equal(test_aia_data, aia_data))
    assert test_aia_header == aia_header


def test_write_file_fits_bytes(tmpdir):
    # Aim is to verify that we can write a FITS file via a file-like object and read it back correctly
    aia_header, aia_data = read_file(TEST_AIA_IMAGE)[0][::-1]
    filepath = tmpdir / "test_aia_171_image_bytes.fits"
    with open(filepath, "wb") as file_obj:
        write_file(file_obj, aia_data, aia_header, filetype='fits')
    assert filepath.exists()
    test_aia_header, test_aia_data = read_file(filepath)[0][::-1]
    assert np.all(np.equal(test_aia_data, aia_data))
    assert test_aia_header == aia_header


def test_missing_file_extension(tmpdir):
    # Aim is read a FITS file without an extension
    aia_header, aia_data = read_file(TEST_AIA_IMAGE)[0][::-1]
    filepath = tmpdir / "test" # No extension
    with open(filepath, 'wb') as file_obj:
        write_file(file_obj, aia_data, aia_header, filetype='fits')
    assert filepath.exists()
    test_aia_header, test_aia_data = read_file(filepath)[0][::-1]
    assert np.all(np.equal(test_aia_data, aia_data))
    assert test_aia_header == aia_header
