
import tempfile

import numpy as np
import pytest

from sunpy.io import ana
from sunpy.tests.helpers import skip_ana

pytestmark = pytest.mark.filterwarnings("ignore::sunpy.util.exceptions.SunpyDeprecationWarning")

img_size = (456, 345)
img_src = np.arange(np.prod(img_size))
img_src.shape = img_size
img_i8 = img_src*2**8/img_src.max()
img_i8 = img_i8.astype(np.int8)
img_i16 = img_src*2**16/img_src.max()
img_i16 = img_i16.astype(np.int16)
img_f32 = img_src*1.0/img_src.max()
img_f32 = img_f32.astype(np.float32)


@skip_ana
def test_i8c():
    # Test int 8 compressed functions
    afilename = tempfile.NamedTemporaryFile().name
    ana.write(afilename, img_i8, 'testcase', 0)
    img_i8c_rec = ana.read(afilename)
    assert np.sum(img_i8c_rec[0][0] - img_i8) == 0


@skip_ana
def test_i8u():
    # Test int 8 uncompressed functions
    afilename = tempfile.NamedTemporaryFile().name
    ana.write(afilename, img_i8, 'testcase', 1)
    img_i8u_rec = ana.read(afilename)
    assert np.sum(img_i8u_rec[0][0] - img_i8) == 0


@skip_ana
def test_i16c():
    # Test int 16 compressed functions
    afilename = tempfile.NamedTemporaryFile().name
    ana.write(afilename, img_i16, 'testcase', 1)
    img_i16c_rec = ana.read(afilename)
    assert np.sum(img_i16c_rec[0][0] - img_i16) == 0


@skip_ana
def test_i16u():
    # Test int 16 uncompressed functions
    afilename = tempfile.NamedTemporaryFile().name
    ana.write(afilename, img_i16, 'testcase', 0)
    img_i16u_rec = ana.read(afilename)
    assert np.sum(img_i16u_rec[0][0] - img_i16) == 0


@skip_ana
def test_f32u():
    # Test float 32 uncompressed functions
    afilename = tempfile.NamedTemporaryFile().name
    ana.write(afilename, img_f32, 'testcase', 0)
    img_f32u_rec = ana.read(afilename)
    assert np.sum(img_f32u_rec[0][0] - img_f32) == 0


@skip_ana
def test_f32c():
    # TODO: Bug with same code. Needs to be tracked down.
    # Test if float 32 compressed functions
    # ana.write('/tmp/pyana-testf32c', img_f32, 1, 'testcase', 0)
    # img_f32c_rec = ana.read('/tmp/pyana-testf32c', 1)
    # assert_(np.sum(img_f32c_rec[0][1]- img_f32) == 0,
    #        msg="Storing 32 bits float data without compression failed (diff: %g)" % (1.0*np.sum(img_f32c_rec[0][1] - img_f32)))
    afilename = tempfile.NamedTemporaryFile().name
    with pytest.raises(RuntimeError):
        ana.write(afilename, img_f32, 'testcase', 1)


@skip_ana
def test_read_memmap():
    # Test to check that passing memmap doesn't raise an error
    afilename = tempfile.NamedTemporaryFile().name
    ana.write(afilename, img_f32, 'testcase', 0)

    # Memmap is not supported by ana
    data_memmap, _ = ana.read(afilename, memmap=True)[0]
    assert data_memmap.base is None

    data, _ = ana.read(afilename, memmap=False)[0]
    assert data.base is None

    assert np.sum(data_memmap - data) == 0
