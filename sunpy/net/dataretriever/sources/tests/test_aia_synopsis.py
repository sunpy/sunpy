import pytest
from hypothesis import given

import astropy.units as u

from sunpy.net import Fido
from sunpy.net import attrs as a
from sunpy.net._attrs import Time
from sunpy.net.dataretriever.client import QueryResponse
from sunpy.net.dataretriever.sources.aia_synopsis import AIASynopsisClient
from sunpy.net.fido_factory import UnifiedResponse
from sunpy.net.tests.strategies import time_attr
from sunpy.time.timerange import TimeRange


@pytest.fixture
def client():
    return AIASynopsisClient()


@pytest.mark.remote_data
@pytest.mark.parametrize(("timerange", "url_start", "url_end"), [
    (Time('2024/4/21', '2024/4/21 00:01'),
     'https://jsoc1.stanford.edu/data/aia/synoptic/2024/04/21/H0000/AIA20240421_0000_0094.fits',
     'https://jsoc1.stanford.edu/data/aia/synoptic/2024/04/21/H0000/AIA20240421_0000_4500.fits'
     ),
    (Time('2024/5/5', '2024/5/5 00:02'),
     'https://jsoc1.stanford.edu/data/aia/synoptic/2024/05/05/H0000/AIA20240505_0000_0094.fits',
     'https://jsoc1.stanford.edu/data/aia/synoptic/2024/05/05/H0000/AIA20240505_0002_1700.fits',
     ),
])
def test_get_url_for_time_range(client, timerange, url_start, url_end):
    qresponse = client.search(timerange)
    urls = [i['url'] for i in qresponse]
    assert isinstance(urls, list)
    assert urls[0] == url_start
    assert urls[-1] == url_end


@given(time_attr())
def test_can_handle_query(time):
    client = AIASynopsisClient()
    # This is to avoid the horrific representation
    ans = client._can_handle_query(time, a.Instrument.aia, a.Level("1.5s"))
    assert ans
    ans = client._can_handle_query(time, a.Instrument.aia)
    assert not ans
    ans = client._can_handle_query(time, a.Instrument.aia, a.Resolution(4.0))
    assert not ans
    ans = client._can_handle_query(time)
    assert not ans


@pytest.mark.remote_data
def test_query(client):
    time_range = TimeRange("2020-01-01", "2020-01-01 00:01")
    query_result = client.search(
        a.Time(time_range.start, time_range.end),
        a.Wavelength(171 * u.angstrom),
        a.Sample(12 * u.hour),
        a.Instrument("AIA"),
        a.Level("1.5s"),
    )
    assert query_result is not None

    qr1 = client.search(Time('2024/8/9', '2024/8/9 00:01'), a.Level("1.5s"))
    assert isinstance(qr1, QueryResponse)
    assert len(qr1) == 9
    assert qr1['Start Time'][0].isot == "2024-08-09T00:00:00.000"
    assert qr1['Start Time'][-1].isot == "2024-08-09T00:00:00.000"


@pytest.mark.remote_data
def test_get(client):
    qr1 = client.search(Time('2024/8/9', '2024/8/9 00:00:30'), a.Level("1.5s"))
    res = client.fetch(qr1)
    assert len(res) == len(qr1)

    res2 = client.fetch(qr1[0])
    assert len(res2) == 1


@pytest.mark.remote_data
def test_fido(client):
    qr = Fido.search(Time('2024/8/9', '2024/8/9 00:00:30'), a.Level("1.5s"))
    client = qr[0].client
    assert isinstance(qr, UnifiedResponse)
    assert isinstance(client, AIASynopsisClient)
    response = Fido.fetch(qr)
    assert len(response) == qr._numfile


def test_attr_reg():
    assert a.Level.onepointfive_s == a.Level('1.5S')


def test_client_repr(client):
    output = str(client)
    assert output[:50] == 'sunpy.net.dataretriever.sources.aia_synopsis.AIASy'
