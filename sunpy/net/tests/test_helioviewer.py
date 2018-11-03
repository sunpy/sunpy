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
        assert type(client.sources['SDO']['AIA']['4500']['sourceId']) is int

    @skip_glymur
    def test_download_jp2(self, client):
        """Tests getJP2Image API method"""
        source_id = client.sources['SOHO']['EIT']['304']['sourceId']
        filepath = client.download_jp2('2020/01/01', source_id)
        map_ = sunpy.map.Map(filepath)
        assert isinstance(map_, sunpy.map.GenericMap)

    @skip_glymur
    def test_download_jp2_directory_not_exist(self, client, tmpdir):
        """Tests getJP2Image API method"""

        source_id = client.sources['SOHO']['EIT']['304']['sourceId']
        filepath = client.download_jp2(
            '2020/01/01',
            source_id,
            directory=os.path.join(str(tmpdir), 'directorynotexist'))

        assert 'directorynotexist' in filepath
