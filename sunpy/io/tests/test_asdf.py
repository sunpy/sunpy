import numpy as np
import pytest

import sunpy.map
from sunpy.data.test import get_test_filepath
from sunpy.io._asdf import get_header, read, write
from sunpy.io._header import FileHeader


@pytest.fixture
def example_map():
    return sunpy.map.Map(get_test_filepath("aia_171_level1.fits"))


@pytest.fixture
def example_map_path(tmp_path, example_map):
    file_path = tmp_path / "map.asdf"
    write(file_path, example_map.data, example_map.meta)
    return file_path


def test_read(tmp_path, example_map, example_map_path):
    cont = read(example_map_path)
    assert isinstance(cont, list)
    assert isinstance(cont[0][0], np.ndarray)
    assert isinstance(cont[0][1], FileHeader)
    assert cont[0][0].shape == example_map.data.shape

def test_write(tmp_path, example_map):
    outfile = tmp_path / "test.asdf"
    write(outfile, example_map.data, example_map.meta)
    assert outfile.exists()
    written_data, written_header = read(outfile)[0]
    assert np.array_equal(example_map.data, written_data)
    assert example_map.meta == written_header
    assert example_map.meta == get_header(outfile)[0]

def test_get_header(example_map_path):
    header = get_header(example_map_path)[0]
    assert isinstance(header, FileHeader)

def test_save(tmpdir, example_map):
    outpath = tmpdir / "fits_map_in_asdf.asdf"
    example_map.save(outpath)
    assert outpath.exists()
    map_obj_asdf = sunpy.map.Map(outpath)
    assert dict(map_obj_asdf.meta) == dict(example_map.meta)
    assert np.array_equal(map_obj_asdf.data, example_map.data)


