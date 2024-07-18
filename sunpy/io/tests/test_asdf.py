import numpy as np

import asdf

import sunpy.map
from sunpy.data.test import get_test_filepath
from sunpy.io._asdf import get_header, get_keys_name, read, write
from sunpy.io._header import FileHeader

map_for_asdf = get_test_filepath("aiamap_genericmap_1.0.0.asdf")

def test_save(tmpdir):
    save_test = sunpy.map.Map(map_for_asdf)
    outfile = tmpdir / 'save_test.asdf'
    save_test.save(str(outfile))
    assert outfile.exists()

def test_read():
    cont = read(map_for_asdf)
    assert isinstance(cont,list)
    assert isinstance(cont[0][0],np.ndarray)
    assert isinstance(cont[0][1],FileHeader)
    assert cont[0][0].shape ==  (2,2)

def test_write(tmpdir):
    data, header = read(map_for_asdf)[0]
    outfile = tmpdir / "test.asdf"
    write(str(outfile), data, header)
    assert outfile.exists()

def test_get_header():
    header = get_header(map_for_asdf)[0]
    assert isinstance(header, FileHeader)
    assert header == sunpy.map.Map(map_for_asdf).meta

def test_keys_name():
    with asdf.open(map_for_asdf) as af:
        rootkeys = af.tree.keys()
        main_data_keys = [key for key in rootkeys if key not in ['asdf_library', 'history']]
        assert get_keys_name(map_for_asdf) == main_data_keys[0]
