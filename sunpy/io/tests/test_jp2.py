"""
JPEG2000 reading tests
"""
from __future__ import absolute_import

#pylint: disable=C0103,R0904,W0201,W0212,W0232,E1103
import numpy as np

from sunpy.data.test import get_test_filepath
from sunpy.io.header import FileHeader
from sunpy.map import GenericMap
from sunpy.map import Map

from sunpy.tests.helpers import skip_glymur

AIA_193_JP2 = get_test_filepath("2013_06_24__17_31_30_84__SDO_AIA_AIA_193.jp2")

@skip_glymur

def test_read_data():
    """Tests the reading of the JP2 data"""
    import glymur
    data = glymur.Jp2k(AIA_193_JP2).read()
    assert isinstance(data, np.ndarray)

@skip_glymur
def test_read_header():
    """Tests the reading of the JP2 header"""
    from sunpy.io.jp2 import get_header
    header = get_header(AIA_193_JP2)[0]
    assert isinstance(header, FileHeader)

@skip_glymur
def test_read_file():
    """Tests the reading of the complete JP2 file and its conversion into a
    SunPy map"""
    map_ = Map(AIA_193_JP2)
    assert isinstance(map_, GenericMap)
