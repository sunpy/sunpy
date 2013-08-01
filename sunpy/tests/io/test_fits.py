from sunpy.io.fits import get_header, extract_waveunit
from sunpy.data.sample import RHESSI_IMAGE, EIT_195_IMAGE, AIA_171_IMAGE,\
    SWAP_LEVEL1_IMAGE


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


def test_extract_waveunit_from_wavelnth_comment():
    # the key WAVEUNIT does not exist, but the comment belonging to the key
    # WAVELNTH can be parsed so that the wave unit can be read
    # SWAP_LEVEL1_IMAGE
    waveunit = extract_waveunit(get_header(SWAP_LEVEL1_IMAGE)[0])
    assert waveunit == 'Angstrom'
