import pytest

import sunpy.net.dataretriever.sources.eve as eve
from sunpy.net import Fido
from sunpy.net import attrs as a
from sunpy.net._attrs import Instrument, Level, Time
from sunpy.net.vso.attrs import Source
from sunpy.net.dataretriever.client import QueryResponse
from sunpy.net.fido_factory import UnifiedResponse
from sunpy.net.vso import VSOClient
from sunpy.time import parse_time
from sunpy.time.timerange import TimeRange

LCClient = eve.EVEClient()


@pytest.mark.remote_data
@pytest.mark.parametrize("timerange,url_start,url_end", [
    (TimeRange('2012/4/21', '2012/4/21'),
     'http://lasp.colorado.edu/eve/data_access/evewebdata/quicklook/L0CS/SpWx/2012/20120421_EVE_L0CS_DIODES_1m.txt',
     'http://lasp.colorado.edu/eve/data_access/evewebdata/quicklook/L0CS/SpWx/2012/20120421_EVE_L0CS_DIODES_1m.txt'
     ),
    (TimeRange('2012/5/5', '2012/5/6'),
     'http://lasp.colorado.edu/eve/data_access/evewebdata/quicklook/L0CS/SpWx/2012/20120505_EVE_L0CS_DIODES_1m.txt',
     'http://lasp.colorado.edu/eve/data_access/evewebdata/quicklook/L0CS/SpWx/2012/20120506_EVE_L0CS_DIODES_1m.txt',
     ),
    (TimeRange('2012/7/7', '2012/7/14'),
     'http://lasp.colorado.edu/eve/data_access/evewebdata/quicklook/L0CS/SpWx/2012/20120707_EVE_L0CS_DIODES_1m.txt',
     'http://lasp.colorado.edu/eve/data_access/evewebdata/quicklook/L0CS/SpWx/2012/20120714_EVE_L0CS_DIODES_1m.txt',
     )
])
def test_get_url_for_time_range(timerange, url_start, url_end):
    urls = LCClient._get_url_for_timerange(timerange)
    assert isinstance(urls, list)
    assert urls[0] == url_start
    assert urls[-1] == url_end


def test_can_handle_query():
    ans1 = eve.EVEClient._can_handle_query(
        Time('2012/8/9', '2012/8/10'), Instrument('eve'), Level(0))
    assert ans1 is True
    ans2 = eve.EVEClient._can_handle_query(Time('2012/7/7', '2012/7/7'))
    assert ans2 is False
    ans3 = eve.EVEClient._can_handle_query(
        Time('2012/8/9', '2012/8/10'), Instrument('eve'), Source('sdo'))
    assert ans3 is False
    ans4 = eve.EVEClient._can_handle_query(
        Time('2012/8/9', '2012/8/10'), Instrument('eve'), Level('0CS'))
    assert ans4 is True
    ans5 = eve.EVEClient._can_handle_query(
        Time('2012/8/9', '2012/8/10'), Instrument('eve'), Level('wibble'))
    assert ans5 is False
    ans6 = eve.EVEClient._can_handle_query(
        Time('2012/8/9', '2012/8/10'), Instrument('eve'), Level(0.5))
    assert ans6 is False


@pytest.mark.remote_data
def test_query():
    qr1 = LCClient.search(Time('2012/8/9', '2012/8/10'), Instrument('eve'))
    assert isinstance(qr1, QueryResponse)
    assert len(qr1) == 2
    assert qr1.time_range().start == parse_time('2012/08/09')
    assert qr1.time_range().end == parse_time('2012/08/11')  # includes end.


@pytest.mark.remote_data
@pytest.mark.parametrize("time,instrument", [
    (Time('2012/11/27', '2012/11/27'), Instrument('eve')),
])
def test_get(time, instrument):
    qr1 = LCClient.search(time, instrument)
    res = LCClient.fetch(qr1)
    assert len(res) == len(qr1)


@pytest.mark.remote_data
@pytest.mark.parametrize(
    'query',
    [(a.Time('2012/10/4', '2012/10/6') & a.Instrument('eve') & a.Level(0))])
def test_fido(query):
    qr = Fido.search(query)
    client = qr.get_response(0).client
    assert isinstance(qr, UnifiedResponse)
    assert isinstance(client, eve.EVEClient)
    response = Fido.fetch(qr)
    assert len(response) == qr._numfile


@pytest.mark.remote_data
@pytest.mark.parametrize(
    'time',
    [(a.Time('2012/10/4', '2012/10/6')), (a.Time('2012/11/27', '2012/11/27'))])
def test_levels(time):
    """
    Test the correct handling of level
    Level 0 comes from EVEClient, other levels from EVE.
    """
    eve_a = a.Instrument('EVE')
    qr = Fido.search(time, eve_a)
    clients = {type(a.client) for a in qr.responses}
    assert clients == {VSOClient}

    qr = Fido.search(time, eve_a, a.Level(0))
    clients = {type(a.client) for a in qr.responses}
    assert clients == {eve.EVEClient}

    # This is broken because the VSO Eve client doesn't provide a way of allowing Level.
    #qr = Fido.search(time, eve_a, a.Level(0) | a.Level(1))
    #clients = {type(a.client) for a in qr.responses}
    #assert clients == {eve.EVEClient}


def test_attr_reg():
    a.Instrument.eve = a.Instrument('EVE')
    a.Level.zero = a.Level('0')
