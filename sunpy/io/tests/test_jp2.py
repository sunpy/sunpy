
import numpy as np

from sunpy.data.test import get_test_filepath
from sunpy.io import _fits, _jp2
from sunpy.io.header import FileHeader
from sunpy.tests.helpers import skip_glymur

AIA_193_JP2 = get_test_filepath("2013_06_24__17_31_30_84__SDO_AIA_AIA_193.jp2")
TEST_AIA_IMAGE = get_test_filepath('aia_171_level1.fits')


@skip_glymur
def test_read():
    data, header = _jp2.read(AIA_193_JP2)[0]
    assert isinstance(data, np.ndarray)
    assert data.shape == (410, 410)
    assert isinstance(header, FileHeader)


@skip_glymur
def test_read_header():
    header = _jp2.get_header(AIA_193_JP2)[0]
    assert isinstance(header, FileHeader)


@skip_glymur
def test_read_memmap():
    data, _ = _jp2.read(AIA_193_JP2, memmap=True)[0]
    # The data is shared, with what, I am unclear
    assert data.base is not None
    # Keyword is not passed in the function call
    data, _ = _jp2.read(AIA_193_JP2, memmap=False)[0]
    assert data.base is not None


@skip_glymur
def test_simple_write(tmpdir):
    data, header = _fits.read(TEST_AIA_IMAGE)[0]
    outfile = tmpdir / "test.jp2"
    _jp2.write(str(outfile), data, header)
    assert outfile.exists()

    # Sanity check that reading back the jp2 returns coherent data
    jp2_readback = _jp2.read(outfile)
    assert header['DATE'] == jp2_readback[0].header['DATE']
