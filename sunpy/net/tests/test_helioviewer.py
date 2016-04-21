"""
Helioviewer Client tests
"""
from __future__ import absolute_import

#pylint: disable=C0103,R0904,W0201,W0212,W0232,E1103
import sunpy
import sunpy.map
import pytest
from sunpy.net.helioviewer import HelioviewerClient

from sunpy.tests.helpers import skip_glymur


@pytest.mark.online
class TestHelioviewerClient:
    """Tests the Helioviewer.org API Client class"""

    def setup_class(self):
        self.client = HelioviewerClient()
        if not self.client.is_online():
            self.__SKIP_TESTS__ = True
            print("Skipping Helioviewer.org tests (server inaccessible)")
        else:
            self.__SKIP_TESTS__ = False

        self.sources = self.client.get_data_sources()

    def teardown_class(self):
        self.client = None

    def test_get_datasources(self):
        """Makes sure datasource query returns a valid result and source id
        is casted to an integer"""
        if self.__SKIP_TESTS__:
            return
        assert type(self.sources['SDO']['AIA']['AIA']['171']['sourceId']) is int

    def test_get_closest_image(self):
        """Tests getClosestImage API method"""
        if self.__SKIP_TESTS__:
            return
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

    @skip_glymur
    def test_download_jp2(self):
        """Tests getJP2Image API method"""
        if self.__SKIP_TESTS__:
            return
        filepath = self.client.download_jp2('2020/01/01', observatory='SOHO',
                                            instrument='MDI', detector='MDI',
                                            measurement='continuum')
        map_ = sunpy.map.Map(filepath)
        assert isinstance(map_, sunpy.map.GenericMap)
