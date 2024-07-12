from sunpy.data.test import get_test_filepath
from sunpy.io._asdf import get_header, read, write
from sunpy.io.header import FileHeader

AIA_MAP = get_test_filepath("aiamap_genericmap_1.0.0.asdf")


def test_read():
    cont = read(AIA_MAP)
    assert isinstance(cont,list)

def test_write(tmpdir):
    data,header = write(AIA_MAP)[0]
    outfile = tmpdir / "test.asdf"
    write(str(outfile),data,header)
    assert outfile.exists()

def test_get_header():
    header = get_header(AIA_MAP)[0]
    assert isinstance(header,FileHeader)
