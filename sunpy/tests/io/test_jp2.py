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

# Tests to see if the main glymur library imports correctly
try:
    __import__("glymur")
except ImportError:
    glymur_imports = True
else:
    glymur_imports = False

@pytest.mark.skipif("glymur_imports is False")
def test_read_data():
    """Tests the reading of the JP2 data"""
    data = glymur.Jp2k(filepath).read()
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
    map_ = Map(filepath)
    assert isinstance(map_, GenericMap)


    @pytest.mark.online
    @pytest.mark.skipif("__SKIP_TESTS__ is True")
    def test_get_closest_image(self):
        """Tests getClosestImage API method"""
        # check basic query
        im1 = self.client.get_closest_image('1994/01/01', 
                                              observatory='SOHO', 
                                              instrument='EIT', 
                                              detector='EIT', 
                                              measurement='195')
        assert im1['width'] == im1['height'] == 1024
        
        # result should be same when using source id to query
        source_id = self.sources['SOHO']['EIT']['EIT']['195']['sourceId']
        
        im2 = self.client.get_closest_image('1994/01/01', sourceId=source_id)
        
        assert im1 == im2
    
    @pytest.mark.online
    @pytest.mark.skipif("__SKIP_TESTS__ is True")
    def test_download_jp2(self):
        """Tests getJP2Image API method"""
        filepath = self.client.download_jp2('2020/01/01', observatory='SOHO', 
                                            instrument='MDI', detector='MDI',
                                            measurement='continuum')
        map_ = sunpy.map.Map(filepath)
        assert isinstance(map_, sunpy.map.GenericMap)
