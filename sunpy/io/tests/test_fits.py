import os
from pathlib import Path
from collections import OrderedDict

import pytest

import astropy.io.fits as fits
from astropy.utils.exceptions import AstropyUserWarning

import sunpy.data.test
import sunpy.io.fits
from sunpy.data.test.waveunit import MEDN_IMAGE, MQ_IMAGE, NA_IMAGE, SVSM_IMAGE
from sunpy.io.fits import extract_waveunit, get_header, header_to_fits
from sunpy.util import MetaDict, SunpyUserWarning

testpath = sunpy.data.test.rootdir

RHESSI_IMAGE = os.path.join(testpath, 'hsi_image_20101016_191218.fits')
EIT_195_IMAGE = os.path.join(testpath, 'EIT/efz20040301.000010_s.fits')
AIA_171_IMAGE = os.path.join(testpath, 'aia_171_level1.fits')
SWAP_LEVEL1_IMAGE = os.path.join(testpath, 'SWAP/resampled1_swap.fits')


def read_hdus():
    pairs = sunpy.io.fits.read(RHESSI_IMAGE)
    assert len(pairs) == 4


def read_hdu_int():
    pairs = sunpy.io.fits.read(RHESSI_IMAGE, hdus=1)
    assert len(pairs) == 1


def read_hdus_list():
    pairs = sunpy.io.fits.read(RHESSI_IMAGE, hdus=[1, 2])
    assert len(pairs) == 2


def read_hdus_gen():
    pairs = sunpy.io.fits.read(RHESSI_IMAGE, hdus=range(0, 1))
    assert len(pairs) == 2


def test_extract_waveunit_missing_waveunit_key_and_missing_wavelnth_comment():
    waveunit = extract_waveunit(get_header(RHESSI_IMAGE)[0])
    assert waveunit is None


def test_missing_waveunit_in_wavelnth_comment():
    # the comment of the key WAVELNTH has the value
    # '171 = Fe IX/X, 195 = Fe XII,' which contains no unit information
    waveunit = extract_waveunit(get_header(EIT_195_IMAGE)[0])
    assert waveunit is None


def test_extract_waveunit_from_waveunit_key():
    # the key WAVEUNIT can be accessed and returned directly
    waveunit = extract_waveunit(get_header(AIA_171_IMAGE)[0])
    assert waveunit == 'angstrom'


def test_extract_waveunit_minus9():
    # value of WAVEUNIT is -9
    with pytest.warns(AstropyUserWarning, match='File may have been truncated'):
        waveunit = extract_waveunit(get_header(MEDN_IMAGE)[0])
    assert waveunit == 'nm'


def test_extract_waveunit_minus10():
    # value of WAVEUNIT is -10
    with pytest.warns(AstropyUserWarning, match='File may have been truncated'):
        waveunit = extract_waveunit(get_header(MQ_IMAGE)[0])
    assert waveunit == 'angstrom'


def test_extract_waveunit_waveunitcomment():
    # comment of WAVEUNIT is: "in meters"
    with pytest.warns(AstropyUserWarning, match='File may have been truncated'):
        waveunit = extract_waveunit(get_header(NA_IMAGE)[0])
    assert waveunit == 'm'


def test_extract_waveunit_wavelnthcomment_brackets():
    # WAVELNTH comment is: "[Angstrom] bandpass peak response"
    waveunit = extract_waveunit(get_header(SWAP_LEVEL1_IMAGE)[0])
    assert waveunit == 'angstrom'


def test_extract_waveunit_wavelnthcomment_parentheses():
    # WAVELNTH comment is: "Observed wavelength (nm)"
    with pytest.warns(AstropyUserWarning, match='File may have been truncated'):
        waveunit = extract_waveunit(get_header(SVSM_IMAGE)[0])
    assert waveunit == 'nm'


def test_simple_write(tmpdir):
    data, header = sunpy.io.fits.read(AIA_171_IMAGE)[0]
    outfile = tmpdir / "test.fits"
    sunpy.io.fits.write(str(outfile), data, header)
    assert outfile.exists()


def test_extra_comment_write(tmpdir):
    data, header = sunpy.io.fits.read(AIA_171_IMAGE)[0]
    header["KEYCOMMENTS"]["TEST"] = "Hello world"
    outfile = tmpdir / "test.fits"
    sunpy.io.fits.write(str(outfile), data, header)
    assert outfile.exists()


def test_simple_write_compressed(tmpdir):
    data, header = sunpy.io.fits.read(AIA_171_IMAGE)[0]
    outfile = tmpdir / "test.fits"
    sunpy.io.fits.write(str(outfile), data, header, hdu_type=fits.CompImageHDU)
    assert outfile.exists()
    with fits.open(str(outfile)) as hdul:
        assert len(hdul) == 2
        assert isinstance(hdul[1], fits.CompImageHDU)


def test_write_with_metadict_header_astropy(tmpdir):
    fits_file = fits.open(AIA_171_IMAGE)
    data, header = fits_file[0].data, fits_file[0].header
    meta_header = MetaDict(OrderedDict(header))
    temp_file = tmpdir / "temp.fits"
    with pytest.warns(SunpyUserWarning, match='The meta key comment is not valid ascii'):
        sunpy.io.fits.write(str(temp_file), data, meta_header)
    assert temp_file.exists()


# Various warnings are thrown in this test, but we just want to check that the code
# works without exceptions
@pytest.mark.filterwarnings('ignore')
def test_fitsheader():
    """Test that all test data can be converted back to a FITS header."""
    extensions = ('fts', 'fits')
    for ext in extensions:
        for ffile in Path(testpath).glob(f"*.{ext}*"):
            fits_file = fits.open(ffile)
            fits_file.verify("fix")
            data, header = fits_file[0].data, fits_file[0].header
            meta_header = MetaDict(OrderedDict(header))
            sunpy.io.fits.header_to_fits(meta_header)
            fits_file.close()


def test_warn_nonascii():
    # Check that a non-ascii character raises a warning and not an error
    with pytest.warns(SunpyUserWarning, match='not valid ascii'):
        fits = header_to_fits({'bad': 'test\t',
                               'good': 'test'})
    assert 'GOOD' in fits.keys()
    assert 'BAD' not in fits.keys()


def test_warn_longkey():
    # Check that a key that is too long raises a warning and not an error
    with pytest.warns(SunpyUserWarning, match='The meta key badlongkey is too long'):
        fits = header_to_fits({'badlongkey': 'test',
                               'goodkey': 'test'})
    assert 'GOODKEY' in fits.keys()
    assert 'BADLONGKEY' not in fits.keys()
