"""
JPEG2000 reading tests
"""
from __future__ import absolute_import

#pylint: disable=C0103,R0904,W0201,W0212,W0232,E1103
import numpy as np
import pytest
from sunpy.data.test import AIA_193_JP2
from sunpy.io.jp2 import get_header
from sunpy.io.header import FileHeader
from sunpy.map import GenericMap
from sunpy.map import Map

# SunPy's JPEG2000 capabilities rely on the glymur library.  First we check to
# make sure that glymur imports correctly before proceeding.
try:
    import glymur
except ImportError:
    glymur_imports = False
else:
    glymur_imports = True

@pytest.mark.skipif("glymur_imports is False")
def test_read_data():
    """Tests the reading of the JP2 data"""
    data = glymur.Jp2k(AIA_193_JP2).read()
    assert isinstance(data, np.ndarray)

@pytest.mark.skipif("glymur_imports is False")
def test_read_header():
    """Tests the reading of the JP2 header"""
    header = get_header(AIA_193_JP2)[0]
    assert isinstance(header, FileHeader)

@pytest.mark.skipif("glymur_imports is False")
def test_read_file():
    """Tests the reading of the complete JP2 file and its conversion into a
    SunPy map"""
    map_ = Map(AIA_193_JP2)
    assert isinstance(map_, GenericMap)
