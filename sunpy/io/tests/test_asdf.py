import numpy as np
import pytest

import asdf

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

def test_read_obj_no_object_key(tmp_path, example_map):
    file_path = tmp_path / "map_no_object_key.asdf"
    af = asdf.AsdfFile({"wrong_key": {"meta": dict(example_map.meta), "data": example_map.data}})
    af.write_to(str(file_path))
    with pytest.raises(KeyError, match="The ASDF does not contain 'object' key"):
        read(file_path)

def test_read_obj_no_meta_data(tmp_path, example_map):
    file_path = tmp_path / "map_no_meta_data.asdf"
    af = asdf.AsdfFile({"object": {"wrong_meta": dict(example_map.meta), "wrong_data": example_map.data}})
    af.write_to(str(file_path))
    with pytest.raises(ValueError, match="The object does not have any meta and data"):
        read(file_path)

def test_read_obj_meta_not_dict(tmp_path, example_map):
    file_path = tmp_path / "map_meta_not_dict.asdf"
    af = asdf.AsdfFile({"object": {"meta": "not_a_dict", "data": example_map.data}})
    af.write_to(str(file_path))
    with pytest.raises(TypeError, match="meta must be a dict not <class 'str'>"):
        read(file_path)

def test_read_obj_data_not_ndarray(tmp_path, example_map):
    file_path = tmp_path / "map_data_not_ndarray.asdf"
    af = asdf.AsdfFile({"object": {"meta": dict(example_map.meta), "data": "not_an_ndarray"}})
    af.write_to(str(file_path))
    with pytest.raises(TypeError, match="data must be a ANDArrayType or numpy ndarray not <class 'str'>"):
        read(file_path)
