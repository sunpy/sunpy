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
def test_search():
    query = [a.solarnet.Dataset.eui_level_2 , a.solarnet.Limit(2) , a.Detector("HRI_EUV")]
    url = Fido.search(*query)
    assert isinstance(url[0],QueryResponseTable)
    assert len(url[0]) == 2


def test_can_handle_query(client):
    assert not client._can_handle_query(a.Time("2020/01/02", "2020/01/03"))
    assert not client._can_handle_query(a.solarnet.Limit(10),a.solarnet.Target.ar)
    assert client._can_handle_query(a.solarnet.Dataset.swap_level_1)


def test_solarnet_attrs(client):
    attrs = client.load_solarnet_values()
    assert a.solarnet.Dataset in attrs.keys()
    assert len(attrs[a.solarnet.Dataset]) > 0


@pytest.mark.remote_data
def test_fetch_return_type():
    qr = Fido.search(a.solarnet.Dataset.swap_level_1 & a.solarnet.Limit(1) & a.solarnet.Tags.moon_transit)
    res = Fido.fetch(qr)
    assert len(qr) == 1
    assert isinstance(res, Results)


@pytest.mark.remote_data
def test_fetch_path_specified(tmpdir):
    query = Fido.search(a.solarnet.Dataset.eui_level_2 , a.solarnet.Limit(2) , a.Detector("HRI_EUV"))
    path = Path(tmpdir)
    files = Fido.fetch(query[0,0] , path=path)
    assert path.exists()
    assert len(files) == 1


@pytest.mark.remote_data
def test_default_limit(client):
    search = client.search(a.solarnet.Dataset.lyra_level_2, a.Wavelength(171*u.AA), a.Time("2020/02/04","2022/02/04"))
    assert len(search) == 20


@pytest.mark.remote_data
def test_complex_query():
    search = Fido.search(a.solarnet.Dataset.lyra_level_2 & a.solarnet.Limit(2) | a.solarnet.Dataset.eui_level_2 & a.solarnet.Limit(3))
    # We have two results as they are not combined
    assert len(search) == 2
    # The first query is limited to 2 results
    assert len(search[0]) == 2
    # The second query is limited to 3 results
    assert len(search[1]) == 3


def test_failed_search_raises(client, monkeypatch):
    def mock_urlopen_fail(url):
        raise OSError("Something went very wrong")

    monkeypatch.setattr("urllib.request.urlopen", mock_urlopen_fail)
    with pytest.raises(OSError, match="Failed to fetch data from"):
        client.search(a.solarnet.Dataset.eui_level_2)
