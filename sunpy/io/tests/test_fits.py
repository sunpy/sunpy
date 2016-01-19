import sunpy.io.fits
from sunpy.io.fits import get_header, extract_waveunit

import sunpy.data.test
import os
testpath = sunpy.data.test.rootdir

RHESSI_IMAGE = os.path.join(testpath, 'hsi_image_20101016_191218.fits')
EIT_195_IMAGE = os.path.join(testpath, 'EIT/efz20040301.000010_s.fits')
AIA_171_IMAGE = os.path.join(testpath, 'aia_171_level1.fits')
SWAP_LEVEL1_IMAGE = os.path.join(testpath, 'SWAP/resampled1_swap.fits')

from sunpy.data.test.waveunit import MEDN_IMAGE, MQ_IMAGE, NA_IMAGE, SVSM_IMAGE

def read_hdus():
    pairs = sunpy.io.fits.read(RHESSI_IMAGE)
    assert len(pairs) == 4

def read_hdu_int():
    pairs = sunpy.io.fits.read(RHESSI_IMAGE, hdus=1)
    assert len(pairs) == 1

def read_hdus_list():
    pairs = sunpy.io.fits.read(RHESSI_IMAGE, hdus=[1,2])
    assert len(pairs) == 2

def read_hdus_gen():
    pairs = sunpy.io.fits.read(RHESSI_IMAGE, hdus=xrange(0,1))
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
    waveunit = extract_waveunit(get_header(MEDN_IMAGE)[0])
    assert waveunit == 'nm'


def test_extract_waveunit_minus10():
    # value of WAVEUNIT is -10
    waveunit = extract_waveunit(get_header(MQ_IMAGE)[0])
    assert waveunit == 'angstrom'


def test_extract_waveunit_waveunitcomment():
    # comment of WAVEUNIT is: "in meters"
    waveunit = extract_waveunit(get_header(NA_IMAGE)[0])
    assert waveunit == 'm'


def test_extract_waveunit_wavelnthcomment_brackets():
    # WAVELNTH comment is: "[Angstrom] bandpass peak response"
    waveunit = extract_waveunit(get_header(SWAP_LEVEL1_IMAGE)[0])
    assert waveunit == 'angstrom'


def test_extract_waveunit_wavelnthcomment_parentheses():
    # WAVELNTH comment is: "Observed wavelength (nm)"
    waveunit = extract_waveunit(get_header(SVSM_IMAGE)[0])
    assert waveunit == 'nm'
