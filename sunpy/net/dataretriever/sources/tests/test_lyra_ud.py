import pytest
from hypothesis import given

import sunpy.net.dataretriever.sources.lyra as lyra
from sunpy.net import Fido
from sunpy.net import attrs as a
from sunpy.net._attrs import Instrument, Time
from sunpy.net.dataretriever.client import QueryResponse
from sunpy.net.fido_factory import UnifiedResponse
from sunpy.net.tests.strategies import range_time
from sunpy.time.timerange import TimeRange


@pytest.fixture
def LCClient():
    return lyra.LYRAClient()


@pytest.mark.remote_data
@pytest.mark.parametrize("timerange,url_start,url_end", [
    (TimeRange('2012/1/7', '2012/1/7'),
     'http://proba2.oma.be/lyra/data/bsd/2012/01/07/lyra_20120107-000000_lev2_std.fits',
     'http://proba2.oma.be/lyra/data/bsd/2012/01/07/lyra_20120107-000000_lev2_std.fits'
     ),
    (TimeRange('2012/12/1', '2012/12/2'),
     'http://proba2.oma.be/lyra/data/bsd/2012/12/01/lyra_20121201-000000_lev2_std.fits',
     'http://proba2.oma.be/lyra/data/bsd/2012/12/02/lyra_20121202-000000_lev2_std.fits'
     ),
    (TimeRange('2012/4/7', '2012/4/14'),
     'http://proba2.oma.be/lyra/data/bsd/2012/04/07/lyra_20120407-000000_lev2_std.fits',
     'http://proba2.oma.be/lyra/data/bsd/2012/04/14/lyra_20120414-000000_lev2_std.fits'
     )
])
def test_get_url_for_time_range(LCClient, timerange, url_start, url_end):
    urls = LCClient._get_url_for_timerange(timerange)
    assert isinstance(urls, list)
    assert urls[0] == url_start
    assert urls[-1] == url_end


@given(range_time('2010-01-06'))
def test_can_handle_query(LCClient, time):
    ans1 = LCClient._can_handle_query(
        time, Instrument('lyra'))
    assert ans1 is True
    ans2 = LCClient._can_handle_query(time)
    assert ans2 is False


@pytest.mark.parametrize("time", [
    Time('2015/8/27', '2015/8/27'),
    Time('2016/2/4', '2016/2/6')])
@pytest.mark.remote_data
def test_query(LCClient, time):
    qr1 = LCClient.search(time, Instrument('lyra'))
    assert isinstance(qr1, QueryResponse)
    assert qr1.time_range().start == time.start
    assert qr1.time_range().end == time.end


@pytest.mark.remote_data
@pytest.mark.parametrize("time,instrument", [
    (Time('2013/8/27', '2013/8/27'), Instrument('lyra'))])
def test_get(LCClient, time, instrument):
    qr1 = LCClient.search(time, instrument)
    download_list = LCClient.fetch(qr1)
    assert len(download_list) == len(qr1)


@pytest.mark.remote_data
@pytest.mark.parametrize(
    "time, instrument",
    [(a.Time('2012/10/4', '2012/10/6'), a.Instrument.lyra)])
def test_fido(time, instrument):
    qr = Fido.search(time, instrument)
    assert isinstance(qr, UnifiedResponse)
    response = Fido.fetch(qr)
    assert len(response) == qr._numfile


def test_attr_reg():
    assert a.Instrument.lyra == a.Instrument('LYRA')
    assert a.Level.one == a.Level('1')
    assert a.Level.two == a.Level('2')
    assert a.Level.three == a.Level('3')


def test_client_repr(LCClient):
    """
    Repr check
    """
    output = str(LCClient)
    assert output[:50] == 'sunpy.net.dataretriever.sources.lyra.LYRAClient\n\nP'
