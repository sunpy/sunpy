import numpy as np
import pytest

from sunpy.io import ana
from sunpy.tests.helpers import skip_ana

pytestmark = [
    skip_ana,
    pytest.mark.filterwarnings("ignore::sunpy.util.exceptions.SunpyDeprecationWarning"),
]


@pytest.fixture(scope="module")
def img_size():
    return (456, 345)


@pytest.fixture(scope="module")
def img_src(img_size):
    img = np.arange(np.prod(img_size))
    img.shape = img_size
    return img


@pytest.fixture(scope="module")
def img_i8(img_src):
    arr = img_src * 2**8 / img_src.max()
    return arr.astype(np.int8)


@pytest.fixture(scope="module")
def img_i16(img_src):
    arr = img_src * 2**16 / img_src.max()
    return arr.astype(np.int16)


@pytest.fixture(scope="module")
def img_f32(img_src):
    arr = img_src * 1.0 / img_src.max()
    return arr.astype(np.float32)


def test_roundtrip_int8_compressed(img_i8, tmp_path):
    p = tmp_path / "i8_compressed.ana"
    ana.write(str(p), img_i8, comments="testcase", compress=True)
    (data, _), = ana.read(str(p))
    np.testing.assert_array_equal(data, img_i8)


def test_roundtrip_int8_uncompressed(img_i8, tmp_path):
    p = tmp_path / "i8_uncompressed.ana"
    ana.write(str(p), img_i8, comments="testcase", compress=False)
    (data, _), = ana.read(str(p))
    np.testing.assert_array_equal(data, img_i8)


def test_roundtrip_int16_compressed(img_i16, tmp_path):
    p = tmp_path / "i16_compressed.ana"
    ana.write(str(p), img_i16, comments="testcase", compress=True)
    (data, _), = ana.read(str(p))
    np.testing.assert_array_equal(data, img_i16)


def test_roundtrip_int16_uncompressed(img_i16, tmp_path):
    p = tmp_path / "i16_uncompressed.ana"
    ana.write(str(p), img_i16, comments="testcase", compress=False)
    (data, _), = ana.read(str(p))
    np.testing.assert_array_equal(data, img_i16)


def test_roundtrip_float32_uncompressed(img_f32, tmp_path):
    p = tmp_path / "f32_uncompressed.ana"
    ana.write(str(p), img_f32, comments="testcase", compress=False)
    (data, _), = ana.read(str(p))
    np.testing.assert_array_equal(data, img_f32)


def test_roundtrip_float32_compressed(img_f32, tmp_path):
    p = tmp_path / "f32_compressed.ana"
    ana.write(str(p), img_f32, comments="testcase", compress=True)
    (data, _), = ana.read(str(p))
    np.testing.assert_array_equal(data, img_f32)


def test_memmap_read_returns_copy(img_f32, tmp_path):
    p = tmp_path / "f32_memmap.ana"
    ana.write(str(p), img_f32, comments="testcase", compress=False)
    # Memmap mode should still return an array that owns its data (no view)
    (data_memmap, _), = ana.read(str(p), memmap=True)
    assert data_memmap.base is None
    (data, _), = ana.read(str(p), memmap=False)
    assert data.base is None
    np.testing.assert_array_equal(data_memmap, data)
