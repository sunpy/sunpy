import numpy as np

import asdf

import sunpy.map
from sunpy.data.test import get_test_filepath
from sunpy.io._asdf import get_header, get_keys_name, read, write
from sunpy.io._header import FileHeader

map_for_asdf = get_test_filepath("aiamap_genericmap_1.0.0.asdf")
fits_map = get_test_filepath("aia_171_level1.fits")

def test_read():
    cont = read(map_for_asdf)
    assert isinstance(cont, list)
    assert isinstance(cont[0][0], np.ndarray)
    assert isinstance(cont[0][1], FileHeader)
    assert cont[0][0].shape == (2, 2)

def test_write_and_verify(tmpdir):
    data, header = read(map_for_asdf)[0]
    outfile = tmpdir / "test.asdf"
    write(str(outfile), data, header)
    assert outfile.exists()
    written_data , written_header = read(str(outfile))[0]
    assert np.array_equal(data,written_data)
    assert header == written_header
    assert header == get_header(str(outfile))[0]

def test_load_save_and_verify(tmpdir):
    data ,header = read(map_for_asdf)[0]
    map_obj = sunpy.map.Map(map_for_asdf)
    outpath = tmpdir / "save.asdf"
    map_obj.save(str(outpath))
    assert outpath.exists()
    assert map_obj.meta == header
    assert np.array_equal(data,map_obj.data)

def verify_data_and_header(tmpdir):
    map_obj_fits = sunpy.map.Map(fits_map)
    outpah = tmpdir/"fits_map_in_asdf.asdf"
    map_obj_fits.save(str(outpah))
    assert outpah.exists()
    map_obj_asdf = sunpy.map.Map(str(outpah))
    assert dict(map_obj_asdf.meta) == dict(map_obj_fits.meta)
    assert np.array_equal(map_obj_asdf.data,map_obj_fits.data)

def test_get_header():
    header = get_header(map_for_asdf)[0]
    assert isinstance(header, FileHeader)

def test_keys_name():
    with asdf.open(map_for_asdf) as af:
        rootkeys = af.tree.keys()
        main_data_keys = [key for key in rootkeys if key not in ['asdf_library', 'history']]
        assert get_keys_name(map_for_asdf) == main_data_keys[0]
