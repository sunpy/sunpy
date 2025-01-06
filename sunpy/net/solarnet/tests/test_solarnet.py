from pathlib import Path

import pytest

import sunpy.net.attrs as a
from sunpy.net.base_client import QueryResponseTable
from sunpy.net.solarnet import SolarnetClient


@pytest.fixture
def client():
    return SolarnetClient()

@pytest.mark.remote_data
def test_search(client):
    query = [a.solarnet.Dataset.eui_level_2 , a.solarnet.Limit(2) , a.Detector("HRI_EUV")]
    url = client.search(*query)
    assert isinstance(url,QueryResponseTable)
    assert len(url) == 2


def test_can_handle_query(client):
    assert not client._can_handle_query(a.Time("2020/01/02", "2020/01/03"))
    assert not client._can_handle_query(a.solarnet.Limit(10))
    assert client._can_handle_query(a.solarnet.Dataset.eui_level_2)

def test_solarnet_attrs(client):
    attrs = client.load_solarnet_values()
    assert a.solarnet.Dataset in attrs.keys()
    assert len(attrs[a.solarnet.Dataset]) > 0

@pytest.mark.remote_data
def test_fetch(client,tmpdir):
    query = client.search(a.solarnet.Dataset.eui_level_2 , a.solarnet.Limit(2) , a.Detector("HRI_EUV"))
    path = Path(tmpdir) / "test_file_1"

    #calling fetch
    client.fetch(query[0],path = str(path))
    assert path.exists()
    expected_file_name = str(query[0]["name"]) + ".fits"
    expected_file = path / expected_file_name
    assert expected_file.exists()
