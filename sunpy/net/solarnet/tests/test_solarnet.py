from pathlib import Path

import pytest
from parfive import Results

import astropy.units as u

import sunpy.net.attrs as a
from sunpy.net import Fido
from sunpy.net.base_client import QueryResponseTable
from sunpy.net.solarnet import SOLARNETClient


@pytest.fixture
def client():
    return SOLARNETClient()


def test_info_url(client):
    assert client.info_url == 'https://solarnet2.oma.be'


@pytest.mark.remote_data
def test_search(client):
    query = [a.solarnet.Dataset.eui_level_2 , a.solarnet.Limit(2) , a.Detector("HRI_EUV")]
    url = client.search(*query)
    assert isinstance(url,QueryResponseTable)
    assert len(url) == 2
    assert "metadata_eui_level_2" in url["datasets"]
    assert "HRI_EUV" in url["detector"]

def test_can_handle_query(client):
    assert not client._can_handle_query(a.Time("2020/01/02", "2020/01/03"))
    assert not client._can_handle_query(a.solarnet.Limit(10))
    assert client._can_handle_query(a.solarnet.Dataset.eui_level_2)


def test_solarnet_attrs(client):
    attrs = client.load_solarnet_values()
    assert a.solarnet.Dataset in attrs.keys()
    assert len(attrs[a.solarnet.Dataset]) > 0

@pytest.mark.remote_data
def test_fetch_return_type():
    qr = Fido.search(a.solarnet.Dataset.eui_level_2 & a.solarnet.Limit(1))
    res = Fido.fetch(qr)
    assert isinstance(res, Results)

@pytest.mark.remote_data
def test_fetch_path_specified(client,tmpdir):
    query = client.search(a.solarnet.Dataset.eui_level_2 , a.solarnet.Limit(2) , a.Detector("HRI_EUV"))
    path = Path(tmpdir) / "test_file_1"

    client.fetch(query[0] , path=path)
    assert path.exists()
    expected_file_name = str(query[0]["name"]) + ".fits"
    expected_file = path / expected_file_name
    assert expected_file.exists()


@pytest.mark.remote_data
def test_default_limit(client):
    search = client.search(a.solarnet.Dataset.lyra_level_2, a.Wavelength(171*u.AA), a.Time("2020/02/04","2022/02/04"))
    assert len(search) == 20


@pytest.mark.remote_data
def test_complex_query():
    search = Fido.search(a.solarnet.Dataset.lyra_level_2 & a.solarnet.Limit(2) | a.solarnet.Dataset.eui_level_2 & a.solarnet.Limit(3))
	# We have two results as they are not combined
    assert len(search) == 2
    # We limited the first query to 2
    assert len(search[0]) == 2
    # We limited the second query to 3
    assert len(search[1]) == 3
    assert "lyra_20100106-000000_lev2_std" in search[0]["name"]
    assert "solo_L2_eui-fsi304-image_20200512T085922556_V06" in search[1]["name"]
    assert "metadata_lyra_level_2" in search[0]["datasets"]
    assert "metadata_eui_level_2"  in search[1]["datasets"]
