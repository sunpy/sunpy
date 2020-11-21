"""
Helioviewer Client tests
"""
import os
import urllib
from collections import OrderedDict

import pytest

import sunpy
import sunpy.map
from sunpy.net.helioviewer import HelioviewerClient
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
        pytest.skip("There was a HTTP error {} {} for "
                    "HelioViewer.".format(e.code, e.args))


@pytest.mark.remote_data
class TestHelioviewerClient:
    """Tests the Helioviewer.org API Client class"""

    def test_get_datasources(self, client):
        """
        Tests get_data_sources and data_sources and that they match.
        """
        assert isinstance(client.data_sources, OrderedDict)
        assert isinstance(client.sources, dict)
        # Rough check that the ordered dict is ordered
        assert list(client.data_sources.values())[0:3] == [0, 1, 2]
        aia_4500_id = client.data_sources['SDO', 'AIA', None, '4500']
        aia_4500_id_copy = client.sources['SDO']['AIA']['4500']['sourceId']
        assert isinstance(aia_4500_id, int)
        assert isinstance(aia_4500_id_copy, int)
        assert aia_4500_id == aia_4500_id_copy

    def test_keyvalue_all(self, client):
        """
        Checks that we raise the correct error for these functions.
        """
        with pytest.raises(KeyError):
            client.get_closest_image("2012/01/01")
        with pytest.raises(KeyError):
            client.download_jp2("2012/01/01")

    def test_get_closest_image(self, client):
        """Tests getClosestImage API method"""
        image_meta = client.get_closest_image('1994/01/01', observatory='SOHO',
                                              instrument='EIT', measurement='304')
        assert isinstance(image_meta, dict)
        assert image_meta['id'] == "1795504"
        assert image_meta['width'] == image_meta['height'] == 1024
        assert image_meta['height'] == image_meta['height'] == 1024
        assert image_meta['name'] == 'EIT 304'
        source_id = client.data_sources['SOHO', 'EIT', None, '304']
        image_meta_id = client.get_closest_image('1994/01/01', source_id=source_id)
        assert image_meta == image_meta_id

    def test_download_jp2(self, client):
        """
        Tests getJP2Image API method.
        """
        filepath = client.download_jp2('2012/01/01', observatory='SOHO',
                                       instrument='MDI', measurement='continuum')
        assert "2011_01_11__22_39_00_000__SOHO_MDI_MDI_continuum.jp2" in filepath
        os.remove(filepath)

    def test_get_jp2_header(self, client):
        """
        Tests getJP2Header API method
        """
        header1 = client.get_jp2_header('1994/01/01', observatory='SOHO',
                                        instrument='EIT', measurement='304')
        header2 = client.get_jp2_header('1994/01/01', jp2_id=1795504)
        assert header1 == header2
        assert len(header1) == len(header2) == 1
        assert ('fits' in header1.keys()) and ('fits' in header2.keys())

    @skip_glymur
    def test_download_jp2_map(self, client):
        """
        Tests getJP2Image API method with Map.
        """
        # TODO: make this a figure test.
        filepath = client.download_jp2('2012/01/01', observatory='SOHO',
                                       instrument='MDI', measurement='continuum')
        map_ = sunpy.map.Map(filepath)
        assert isinstance(map_, sunpy.map.GenericMap)
        os.remove(filepath)

    def test_download_directory_not_exist_all(self, client, tmpdir):
        """
        Tests for missing directory.
        """
        fake_dir = os.path.join(str(tmpdir), 'directorynotexist')
        filepath = client.download_jp2('2020/01/01', observatory='SOHO',
                                       instrument='MDI', measurement='continuum',
                                       directory=fake_dir)
        assert 'directorynotexist' in filepath
        os.remove(filepath)
        fake_dir = os.path.join(str(tmpdir), 'directorynotexist_2')
        filepath = client.download_png('2020/01/01', 2.4, "[SOHO,MDI,continuum,1,100]",
                                       directory=fake_dir)
        assert 'directorynotexist_2' in filepath
        os.remove(filepath)

    def test_overwrite_jp2(self, client):
        """
        Tests for that overwrites, overwrites jp2 edition.
        """
        filepath = client.download_jp2('2020/01/01', observatory='SOHO',
                                       instrument='MDI', measurement='continuum',
                                       overwrite=False)
        filepath_2 = client.download_jp2('2020/01/01', observatory='SOHO',
                                         instrument='MDI', measurement='continuum',
                                         overwrite=False)
        assert filepath_2 == filepath
        filepath_3 = client.download_jp2('2020/01/01', observatory='SOHO',
                                         instrument='MDI', measurement='continuum',
                                         overwrite=True)
        assert filepath_3 == filepath

    def test_overwrite_png(self, client):
        """
        Tests for that overwrites, overwrites png edition.
        """
        filepath = client.download_png('2020/01/01', 2.4, "[SOHO,MDI,continuum,1,100]",
                                       overwrite=False)
        filepath_2 = client.download_png('2020/01/01', 2.4, "[SOHO,MDI,continuum,1,100]",
                                         overwrite=False)
        assert filepath_2 is not filepath
        filepath_3 = client.download_png('2020/01/01', 2.4, "[SOHO,MDI,continuum,1,100]",
                                         overwrite=True)
        assert filepath_3 == filepath
