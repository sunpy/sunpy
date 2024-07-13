import numpy as np

import asdf

from sunpy.data.test import get_test_filepath
from sunpy.io._asdf import get_header, get_keys_name, read, write
from sunpy.io._header import FileHeader

AIA_MAP = get_test_filepath("aiamap_genericmap_1.0.0.asdf")


def test_read():
    cont = read(AIA_MAP)
    assert isinstance(cont,list)
    assert isinstance(cont[0][0],np.ndarray)
    assert isinstance(cont[0][1],FileHeader)
    assert cont[0][0].shape ==  (2,2)

def test_write(tmpdir):
    data,header = read(AIA_MAP)[0]
    outfile = tmpdir / "test.asdf"
    write(str(outfile),data,header)
    assert outfile.exists()

def test_get_header():
    header = get_header(AIA_MAP)[0]
    assert isinstance(header,FileHeader)

def test_keys_name():
    with asdf.open(AIA_MAP) as af:
        rootkeys = af.tree.keys()
        main_data_keys = [key for key in rootkeys if key not in ['asdf_library', 'history']]
        assert get_keys_name(AIA_MAP) == main_data_keys[0]
