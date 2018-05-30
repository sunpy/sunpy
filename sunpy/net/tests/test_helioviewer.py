"""
Helioviewer Client tests
"""
from __future__ import absolute_import

import os

import sunpy
import sunpy.map
import pytest
from sunpy.net.helioviewer import HelioviewerClient
from sunpy.extern.six.moves import urllib

from sunpy.tests.helpers import skip_glymur


@pytest.fixture(scope="function")
def client():
    """
    Fixture to create a client and skip tests if not available
    """
    try:
        client = HelioviewerClient()
        client.sources = client.get_data_sources()
        return client
    except urllib.error.HTTPError as e:
        print("There's a HTTP problem {} {}".format(e.code, e.args))
        pytest.skip("HTTP error {}".format(e.code))


@pytest.mark.remote_data
class TestHelioviewerClient:
    """Tests the Helioviewer.org API Client class"""

    def test_get_datasources(self, client):
        """Makes sure datasource query returns a valid result and source id
        is casted to an integer"""
        assert type(client.sources['SDO']['AIA']['AIA']['171']['sourceId']) is int

    def test_get_closest_image(self, client):
        """Tests getClosestImage API method"""
        # check basic query
        im1 = client.get_closest_image('1994/01/01',
                                       observatory='SOHO',
                                       instrument='EIT',
                                       detector='EIT',
                                       measurement='195')
        assert im1['width'] == im1['height'] == 1024

        # result should be same when using source id to query
        source_id = client.sources['SOHO']['EIT']['EIT']['195']['sourceId']

        im2 = client.get_closest_image('1994/01/01', sourceId=source_id)

        assert im1 == im2

    @skip_glymur
    def test_download_jp2(self, client):
        """Tests getJP2Image API method"""
        filepath = client.download_jp2('2020/01/01', observatory='SOHO',
                                       instrument='MDI', detector='MDI',
                                       measurement='continuum')
        map_ = sunpy.map.Map(filepath)
        assert isinstance(map_, sunpy.map.GenericMap)

    @skip_glymur
    def test_download_jp2_directory_not_exist(self, client, tmpdir):
        """Tests getJP2Image API method"""

        filepath = client.download_jp2(
            '2020/01/01',
            observatory='SOHO',
            instrument='MDI',
            detector='MDI',
            measurement='continuum',
            directory=os.path.join(str(tmpdir), 'directorynotexist'))

        assert 'directorynotexist' in filepath
