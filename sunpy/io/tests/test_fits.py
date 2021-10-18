import os
from pathlib import Path
from collections import OrderedDict

import numpy as np
import pytest

import astropy.io.fits as fits
from astropy.utils.exceptions import AstropyUserWarning

import sunpy.data.test
import sunpy.io.fits
from sunpy.data.test.waveunit import MEDN_IMAGE, MQ_IMAGE, NA_IMAGE, SVSM_IMAGE
from sunpy.io.fits import extract_waveunit, get_header, header_to_fits
from sunpy.util import MetaDict, SunpyMetadataWarning

testpath = sunpy.data.test.rootdir

RHESSI_IMAGE = os.path.join(testpath, 'hsi_image_20101016_191218.fits')
EIT_195_IMAGE = os.path.join(testpath, 'EIT/efz20040301.000010_s.fits')
AIA_171_IMAGE = os.path.join(testpath, 'aia_171_level1.fits')
SWAP_LEVEL1_IMAGE = os.path.join(testpath, 'SWAP/resampled1_swap.fits')


# Some of the tests iamges contain an invalid BLANK keyword; ignore the warning
# raised by this
pytestmark = pytest.mark.filterwarnings("ignore:Invalid 'BLANK' keyword in header")


@pytest.mark.parametrize(
    'fname, hdus, length',
    [(RHESSI_IMAGE, None, 4),
     (RHESSI_IMAGE, 1, 1),
     (RHESSI_IMAGE, [1, 2], 2),
     (RHESSI_IMAGE, range(0, 1), 2)]
)
def read_hdus(fname, hdus, length):
    pairs = sunpy.io.fits.read(fname, hdus=hdus)
    assert len(pairs) == length


@pytest.mark.parametrize(
    'fname, waveunit, warn',
    [(RHESSI_IMAGE, None, False),
     (EIT_195_IMAGE, None, False),
     (AIA_171_IMAGE, 'angstrom', False),
     (MEDN_IMAGE, 'nm', True),
     (MQ_IMAGE, 'angstrom', True),
     (NA_IMAGE, 'm', True),
     (SWAP_LEVEL1_IMAGE, 'angstrom', False),
     (SVSM_IMAGE, 'nm', True)]
)
def test_extract_waveunit(fname, waveunit, warn):
    if warn:
        with pytest.warns(AstropyUserWarning, match='File may have been truncated'):
            waveunit = extract_waveunit(get_header(fname)[0])
    else:
        waveunit = extract_waveunit(get_header(fname)[0])
    assert waveunit is waveunit


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


def test_simple_write_compressed_difftypeinst(tmpdir):
    # `hdu_type=fits.CompImageHDU` and `hdu_type=fits.CompImageHDU()`
    # should produce identical FITS files
    data, header = sunpy.io.fits.read(AIA_171_IMAGE)[0]
    outfile_type = str(tmpdir / "test_type.fits")
    outfile_inst = str(tmpdir / "test_inst.fits")
    sunpy.io.fits.write(outfile_type, data, header, hdu_type=fits.CompImageHDU)
    sunpy.io.fits.write(outfile_inst, data, header, hdu_type=fits.CompImageHDU())
    assert fits.FITSDiff(outfile_type, outfile_inst, ignore_comments=['PCOUNT']).identical


@pytest.mark.parametrize(
    'kwargs, should_fail',
    [({}, False),
     ({'quantize_level': -32}, True)]
)
def test_simple_write_compressed_instance(tmpdir, kwargs, should_fail):
    data, header = sunpy.io.fits.read(AIA_171_IMAGE)[0]
    outfile = tmpdir / "test.fits"

    # Ensure HDU instance is used correctly
    hdu = fits.CompImageHDU(data=np.array([0.]), **kwargs)
    hdu.header['HELLO'] = 'world'  # should be in the written file
    hdu.header['TELESCOP'] = 'other'  # should be replaced with 'SDO/AIA'
    hdu.header['NAXIS'] = 5  # should be replaced with 2
    sunpy.io.fits.write(str(outfile), data, header, hdu_type=hdu)
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
    with fits.open(AIA_171_IMAGE) as fits_file:
        data, header = fits_file[0].data, fits_file[0].header
    meta_header = MetaDict(OrderedDict(header))
    temp_file = tmpdir / "temp.fits"
    with pytest.warns(SunpyMetadataWarning, match='The meta key comment is not valid ascii'):
        sunpy.io.fits.write(str(temp_file), data, meta_header)
    assert temp_file.exists()
    fits_file.close()

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
            meta_header = MetaDict(OrderedDict(fits_file[0].header))
            sunpy.io.fits.header_to_fits(meta_header)
            fits_file.close()


def test_warn_nonascii():
    # Check that a non-ascii character raises a warning and not an error
    with pytest.warns(SunpyMetadataWarning, match='not valid ascii'):
        fits = header_to_fits({'bad': 'test\t',
                               'good': 'test'})
    assert 'GOOD' in fits.keys()
    assert 'BAD' not in fits.keys()


def test_warn_nan():
    # Check that a NaN value raises a warning and not an error
    with pytest.warns(SunpyMetadataWarning, match='has a NaN value'):
        fits = header_to_fits({'bad': float('nan'),
                               'good': 1.0})
    assert 'GOOD' in fits.keys()
    assert 'BAD' not in fits.keys()


def test_warn_longkey():
    # Check that a key that is too long raises a warning and not an error
    with pytest.warns(SunpyMetadataWarning, match='The meta key badlongkey is too long'):
        fits = header_to_fits({'badlongkey': 'test',
                               'goodkey': 'test'})
    assert 'GOODKEY' in fits.keys()
    assert 'BADLONGKEY' not in fits.keys()
