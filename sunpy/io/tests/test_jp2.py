import numpy as np

from sunpy.data.test import get_test_filepath
from sunpy.io.header import FileHeader
from sunpy.io.jp2 import get_header, read
from sunpy.tests.helpers import skip_glymur

AIA_193_JP2 = get_test_filepath("2013_06_24__17_31_30_84__SDO_AIA_AIA_193.jp2")


@skip_glymur
def test_read():
    data, header = read(AIA_193_JP2)[0]
    assert isinstance(data, np.ndarray)
    assert data.shape == (4096, 4096)
    assert isinstance(header, FileHeader)


@skip_glymur
def test_read_header():
    header = get_header(AIA_193_JP2)[0]
    assert isinstance(header, FileHeader)


@skip_glymur
def test_read_memmap():
    data, _ = read(AIA_193_JP2, memmap=True)[0]
    # The data is shared, with what, I am unclear
    assert data.base is not None
    # Keyword is not passed in the function call
    data, _ = read(AIA_193_JP2, memmap=False)[0]
    assert data.base is not None
