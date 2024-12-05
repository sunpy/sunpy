import mmap
from pathlib import Path
from collections import OrderedDict

import numpy as np
import pytest

import astropy.io.fits as fits

from sunpy.data.test import get_test_data_filenames, get_test_filepath
from sunpy.data.test.waveunit import MEDN_IMAGE, MQ_IMAGE, NA_IMAGE, SVSM_IMAGE
from sunpy.io import _fits
from sunpy.util import MetaDict, SunpyMetadataWarning

TEST_RHESSI_IMAGE = get_test_filepath('hsi_image_20101016_191218.fits')
TEST_AIA_IMAGE = get_test_filepath('aia_171_level1.fits')
TEST_EIT_HEADER = get_test_filepath('EIT_header/efz20040301.000010_s.header')
TEST_SWAP_HEADER = get_test_filepath('SWAP/resampled1_swap.header')
TEST_GONG_HEADER = get_test_filepath('gong_synoptic.header')
# Some of the tests images contain an invalid BLANK keyword
# ignore the warning raised by this
pytestmark = pytest.mark.filterwarnings("ignore:Invalid 'BLANK' keyword in header")


@pytest.mark.parametrize(
    ('fname', 'hdus', 'length'),
    [(TEST_RHESSI_IMAGE, None, 4),
     (TEST_RHESSI_IMAGE, 1, 1),
     (TEST_RHESSI_IMAGE, [1, 2], 2),
     (TEST_RHESSI_IMAGE, range(0, 2), 2)]
)
def test_read_hdus(fname, hdus, length):
    pairs = _fits.read(fname, hdus=hdus)
    assert len(pairs) == length


@pytest.mark.parametrize(
    ('fname', 'waveunit'),
    [(TEST_RHESSI_IMAGE, None),
     (TEST_EIT_HEADER, None),
     (TEST_AIA_IMAGE, 'angstrom'),
     (MEDN_IMAGE, 'nm'),
     (MQ_IMAGE, 'angstrom'),
     (NA_IMAGE, 'm'),
     (TEST_SWAP_HEADER, 'angstrom'),
     (SVSM_IMAGE, 'nm'),
     (TEST_GONG_HEADER, 'nm'),]
)
def test_extract_waveunit(fname, waveunit):
    if Path(fname).suffix == '.header':
        header = _fits.format_comments_and_history(fits.Header.fromtextfile(fname))
    else:
        header = _fits.get_header(fname)[0]
    waveunit = _fits.extract_waveunit(header)
    assert waveunit is waveunit


def test_simple_write(tmpdir):
    data, header = _fits.read(TEST_AIA_IMAGE)[0]
    outfile = tmpdir / "test.fits"
    _fits.write(str(outfile), data, header)
    assert outfile.exists()


def test_extra_comment_write(tmpdir):
    data, header = _fits.read(TEST_AIA_IMAGE)[0]
    header["KEYCOMMENTS"]["TEST"] = "Hello world"
    outfile = tmpdir / "test.fits"
    _fits.write(str(outfile), data, header)
    assert outfile.exists()


def test_simple_write_compressed(tmpdir):
    data, header = _fits.read(TEST_AIA_IMAGE)[0]
    outfile = tmpdir / "test.fits"
    _fits.write(str(outfile), data, header, hdu_type=fits.CompImageHDU)
    assert outfile.exists()
    with fits.open(str(outfile)) as hdul:
        assert len(hdul) == 2
        assert isinstance(hdul[1], fits.CompImageHDU)


def test_simple_write_compressed_difftypeinst(tmpdir):
    # `hdu_type=fits.CompImageHDU` and `hdu_type=fits.CompImageHDU()`
    # should produce identical FITS files
    data, header = _fits.read(TEST_AIA_IMAGE)[0]
    outfile_type = str(tmpdir / "test_type.fits")
    outfile_inst = str(tmpdir / "test_inst.fits")
    _fits.write(outfile_type, data, header, hdu_type=fits.CompImageHDU, output_verify="silentfix")
    _fits.write(outfile_inst, data, header, hdu_type=fits.CompImageHDU(), output_verify="silentfix")
    assert fits.FITSDiff(outfile_type, outfile_inst, ignore_comments=["SIMPLE"]).identical


@pytest.mark.parametrize(
    ('kwargs', 'should_fail'),
    [({}, False),
     ({'quantize_level': -32}, True)]
)
def test_simple_write_compressed_instance(tmpdir, kwargs, should_fail):
    data, header = _fits.read(TEST_AIA_IMAGE)[0]
    outfile = tmpdir / "test.fits"

    # Ensure HDU instance is used correctly
    hdu = fits.CompImageHDU(data=np.array([0.]), **kwargs)
    hdu.header['HELLO'] = 'world'  # should be in the written file
    hdu.header['TELESCOP'] = 'other'  # should be replaced with 'SDO/AIA'
    hdu.header['NAXIS'] = 5  # should be replaced with 2
    _fits.write(str(outfile), data, header, hdu_type=hdu, output_verify="silentfix")
    assert outfile.exists()
    with fits.open(str(outfile)) as hdul:
        assert len(hdul) == 2
        assert isinstance(hdul[1], fits.CompImageHDU)
        assert hdul[1].header['HELLO'] == 'world'
        assert hdul[1].header['TELESCOP'] == 'SDO/AIA'
        assert hdul[1].header['NAXIS'] == 2
        data_preserved = hdul[1].data == pytest.approx(data, abs=10)
        print(np.abs(hdul[1].data - data).max())
        print(kwargs)
        if should_fail:  # high compression setting preserved
            assert not data_preserved
        else:
            assert data_preserved


def test_write_with_metadict_header_astropy(tmpdir):
    with fits.open(TEST_AIA_IMAGE) as fits_file:
        data, header = fits_file[0].data, fits_file[0].header
    meta_header = MetaDict(OrderedDict(header))
    temp_file = tmpdir / "temp.fits"
    _fits.write(str(temp_file), data, meta_header)
    assert temp_file.exists()
    fits_file.close()

# Various warnings are thrown in this test, but we just want to check that the code
# works without exceptions


@pytest.mark.filterwarnings('ignore')
def test_fitsheader():
    """Test that all test data can be converted back to a FITS header."""
    extensions = ('.fts', '.fits')
    for ext in extensions:
        bad_file = "not_actually_fits.fits"
        test_files = [f for f in get_test_data_filenames() if f.suffix == ext and f.name != bad_file]
        for ffile in test_files:
            fits_file = fits.open(ffile)
            fits_file.verify("fix")
            meta_header = MetaDict(OrderedDict(fits_file[0].header))
            _fits.header_to_fits(meta_header)
            fits_file.close()


def test_warn_nonascii():
    # Check that a non-ascii character raises a warning and not an error
    with pytest.warns(SunpyMetadataWarning, match='not valid ascii'):
        fits = _fits.header_to_fits({'bad': 'test\t',
                                     'good': 'test'})
    assert 'GOOD' in fits.keys()
    assert 'BAD' not in fits.keys()


def test_warn_nan():
    # Check that a NaN value raises a warning and not an error
    with pytest.warns(SunpyMetadataWarning, match='has a NaN value'):
        fits = _fits.header_to_fits({'bad': float('nan'),
                                     'good': 1.0})
    assert 'GOOD' in fits.keys()
    assert 'BAD' not in fits.keys()


def test_warn_longkey():
    # Check that a key that is too long raises a warning and not an error
    with pytest.warns(SunpyMetadataWarning, match='The meta key badlongkey is too long'):
        fits = _fits.header_to_fits({'badlongkey': 'test',
                                     'goodkey': 'test'})
    assert 'GOODKEY' in fits.keys()
    assert 'BADLONGKEY' not in fits.keys()


def test_read_memmap():
    data, _ = _fits.read(TEST_AIA_IMAGE, memmap=True)[0]
    assert data.base is not None
    assert isinstance(data.base, mmap.mmap)

    data, _ = _fits.read(TEST_AIA_IMAGE, memmap=False)[0]
    assert data.base is None


def test_merge_multiple_comment_history_entries():
    header = fits.Header()
    history = "First history entry\nSecond history entry"
    comment = "First comment\nSecond comment"
    for hist in history.split("\n"):
        header.add_history(hist)
    for com in comment.split("\n"):
        header.add_comment(com)
    file_header = _fits.format_comments_and_history(header)
    assert file_header["HISTORY"] == history
    assert file_header["COMMENT"] == comment


def test_split_multiline_comment_history_entries():
    header_dict = {
        "history": "First history entry\nSecond history entry",
        "comment": "First comment\nSecond comment",
    }
    fits_header = _fits.header_to_fits(header_dict)
    for k, v in header_dict.items():
        for v_fits, v_dict in zip(fits_header[k.upper()], v.split("\n")):
            assert v_fits == v_dict


def test_roundtrip_comment_history_entries(tmpdir):
    # Construct FITS file from header with multiple history and comment cards
    header = fits.Header.fromtextfile(get_test_filepath("solo_L1_eui-fsi304-image_20201021T145510206_V03.header"))
    data = np.random.rand(header['NAXIS2'], header['NAXIS1'])
    outfilename = str(tmpdir / 'roundtrip-comments-history-astropy.fits')
    fits.writeto(outfilename, data, header)
    # Check roundtripping
    (dh_pair, ) = _fits.read(outfilename)
    header_dict = dh_pair.header
    assert header_dict["HISTORY"] == "\n".join(header["HISTORY"])
    assert header_dict["COMMENT"] == "\n".join(header["COMMENT"])
    outfilename = str(tmpdir / 'roundtrip-comments-history-sunpy.fits')
    _fits.write(outfilename, data, header_dict)
    with fits.open(outfilename) as hdul:
        header_rt = hdul[0].header
    for key in ("HISTORY", "COMMENT"):
        assert len(header_rt[key]) == len(header[key])
        for card_rt, card in zip(header_rt[key], header[key]):
            assert card_rt == card
