import os
from pathlib import Path

import pytest

import astropy.units as u

import sunpy.net.attrs as a
from sunpy.net import Fido
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
    client.fetch(query[0],path =path)
    assert path.exists()
    expected_file_name = str(query[0]["name"]) + ".fits"
    expected_file = path / expected_file_name
    assert expected_file.exists()
    os.remove(expected_file)

    # Verify the file has been deleted
    assert not expected_file.exists()

@pytest.mark.remote_data
def test_default_limit(client):
    query = [a.solarnet.Dataset.lyra_level_2, a.Wavelength(171*u.AA),a.Time("2020/02/04","2022/02/04")]
    url = client.search(*query)
    assert len(url) == 20

@pytest.mark.remote_data
def test_complex_query():
    query = [a.solarnet.Dataset.lyra_level_2 & a.solarnet.Limit(2) | a.solarnet.Dataset.eui_level_2 & a.solarnet.Limit(3)]
    url = Fido.search(*query)
    assert len(url) == 2
    assert len(url[0]) == 2
    assert len(url[1]) == 3
