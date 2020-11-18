import pytest

import sunpy.net.dataretriever.sources.gong as gong
from sunpy.net import Fido
from sunpy.net import attrs as a
from sunpy.net._attrs import Instrument, Time
from sunpy.net.dataretriever.client import QueryResponse
from sunpy.net.fido_factory import UnifiedResponse


@pytest.fixture
def GSClient():
    return gong.GONGClient()


@pytest.mark.remote_data
@pytest.mark.parametrize("timerange,url_start,url_end", [
    (a.Time('2020/1/30', '2020/2/1'),
     'https://gong2.nso.edu/oQR/zqs/202001/mrzqs200130/mrzqs200130t0004c2227_349.fits.gz',
     'https://gong2.nso.edu/oQR/zqs/202001/mrzqs200131/mrzqs200131t2314c2227_323.fits.gz'
     ),
    (a.Time('2020/4/21', '2020/4/22'),
     'https://gong2.nso.edu/oQR/zqs/202004/mrzqs200421/mrzqs200421t0004c2230_348.fits.gz',
     'https://gong2.nso.edu/oQR/zqs/202004/mrzqs200421/mrzqs200421t2314c2230_335.fits.gz'
     ),
    (a.Time('2006/9/19', '2006/9/19 22:00'),
     'https://gong2.nso.edu/oQR/zqs/200609/mrzqs060919/mrzqs060919t1154c2048_323.fits.gz',
     'https://gong2.nso.edu/oQR/zqs/200609/mrzqs060919/mrzqs060919t1754c2048_320.fits.gz'
     )
])
def test_get_url_for_time_range(GSClient, timerange, url_start, url_end):
    qresponse = GSClient.search(timerange)
    urls = [i['url'] for i in qresponse]
    assert isinstance(urls, list)
    assert urls[0] == url_start
    assert urls[-1] == url_end


@pytest.mark.parametrize(
    "query, result",
    [((a.Time('2020/1/1', '2020/1/2'), a.Instrument('gong')), True),
     ((a.Time('2020/1/1', '2020/1/2'), a.Instrument('goes')), False),
     ((a.Time('2020/1/1', '2020/1/2'), a.Instrument('gong'), a.Physobs('LOS_MAGNETIC_FIELD'),
       a.ExtentType('synoptic')), True),
     ((a.Time('2020/1/1', '2020/1/2'), a.Instrument('gong'),
       a.ExtentType('synoptic')), True),
     ((a.Time('2020/1/1', '2020/1/2'), a.Instrument('gong'),
       a.ExtentType('FULL_DISK')), False)])
def test_can_handle_query(query, result):
    assert gong.GONGClient._can_handle_query(*query) == result


@pytest.mark.remote_data
@pytest.mark.parametrize("time,instrument", [
    (Time('2013/8/27', '2013/8/27'), Instrument('gong')),
    (Time('2020/4/23 17:00', '2020/4/23 21:00'), Instrument('gong')),
])
def test_get(GSClient, time, instrument):
    qr1 = GSClient.search(time, instrument)
    assert isinstance(qr1, QueryResponse)
    download_list = GSClient.fetch(qr1)
    assert len(download_list) == len(qr1)


@pytest.mark.remote_data
@pytest.mark.parametrize(
    "time, instrument",
    [(a.Time('2019/10/4', '2019/10/4 2:00'), a.Instrument('gong')),
     (a.Time('2019/12/31 21:00', '2020/1/1'), a.Instrument('gong'))])
def test_fido(time, instrument):
    qr = Fido.search(time, instrument)
    assert isinstance(qr, UnifiedResponse)
    response = Fido.fetch(qr)
    assert len(response) == qr._numfile
    assert all(map(lambda x: x.endswith('.gz'), response))


def test_attr_reg():
    assert a.Instrument.gong == a.Instrument("GONG")
    assert a.ExtentType.synoptic == a.ExtentType("SYNOPTIC")


def test_client_repr(GSClient):
    """
    Repr check
    """
    output = str(GSClient)
    assert output[:50] == 'sunpy.net.dataretriever.sources.gong.GONGClient\n\nP'
