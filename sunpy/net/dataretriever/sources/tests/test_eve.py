import pytest

import sunpy.net.dataretriever.sources.eve as eve
from sunpy.net import Fido
from sunpy.net import attrs as a
from sunpy.net._attrs import Instrument, Level, Time
from sunpy.net.dataretriever.client import QueryResponse
from sunpy.net.fido_factory import UnifiedResponse
from sunpy.net.vso import VSOClient
from sunpy.time import parse_time


@pytest.fixture
def LCClient():
    return eve.EVEClient()


@pytest.mark.remote_data
@pytest.mark.parametrize(("timerange", "url_start", "url_end"), [
    (Time('2012/4/21', '2012/4/21'),
     'https://lasp.colorado.edu/eve/data_access/eve_data/quicklook/L0CS/SpWx/2012/20120421_EVE_L0CS_DIODES_1m.txt',
     'https://lasp.colorado.edu/eve/data_access/eve_data/quicklook/L0CS/SpWx/2012/20120421_EVE_L0CS_DIODES_1m.txt'
     ),
    (Time('2012/5/5', '2012/5/6'),
     'https://lasp.colorado.edu/eve/data_access/eve_data/quicklook/L0CS/SpWx/2012/20120505_EVE_L0CS_DIODES_1m.txt',
     'https://lasp.colorado.edu/eve/data_access/eve_data/quicklook/L0CS/SpWx/2012/20120506_EVE_L0CS_DIODES_1m.txt',
     ),
    (Time('2012/7/7', '2012/7/14'),
     'https://lasp.colorado.edu/eve/data_access/eve_data/quicklook/L0CS/SpWx/2012/20120707_EVE_L0CS_DIODES_1m.txt',
     'https://lasp.colorado.edu/eve/data_access/eve_data/quicklook/L0CS/SpWx/2012/20120714_EVE_L0CS_DIODES_1m.txt',
     )
])
def test_get_url_for_time_range(LCClient, timerange, url_start, url_end):
    qresponse = LCClient.search(timerange)
    urls = [i['url'] for i in qresponse]
    assert isinstance(urls, list)
    assert urls[0] == url_start
    assert urls[-1] == url_end


def test_can_handle_query(LCClient):
    ans1 = LCClient._can_handle_query(
        Time('2012/8/9', '2012/8/10'), Instrument('eve'), Level(0))
    assert ans1 is True
    ans2 = LCClient._can_handle_query(Time('2012/7/7', '2012/7/7'))
    assert ans2 is False
    ans3 = LCClient._can_handle_query(
        Time('2012/8/9', '2012/8/10'), Instrument('eve'), a.Source('sdo'))
    assert ans3 is False
    ans4 = LCClient._can_handle_query(
        Time('2012/8/9', '2012/8/10'), Instrument('eve'), Level('0CS'))
    assert ans4 is False
    ans5 = LCClient._can_handle_query(
        Time('2012/8/9', '2012/8/10'), Instrument('eve'), Level('wibble'))
    assert ans5 is False
    ans6 = LCClient._can_handle_query(
        Time('2012/8/9', '2012/8/10'), Instrument('eve'), Level(0.5))
    assert ans6 is False


@pytest.mark.remote_data
def test_query(LCClient):
    qr1 = LCClient.search(Time('2012/8/9', '2012/8/10'), Instrument('eve'))
    assert isinstance(qr1, QueryResponse)
    assert len(qr1) == 2
    assert qr1['Start Time'][0].datetime == parse_time('2012/08/09').datetime
    assert qr1['End Time'][1].datetime == parse_time('2012/08/10 23:59:59.999').datetime


@pytest.mark.remote_data
@pytest.mark.parametrize(("time", "instrument"), [
    (Time('2012/11/27', '2012/11/27'), Instrument('eve')),
])
def test_get(LCClient, time, instrument):
    qr1 = LCClient.search(time, instrument)
    res = LCClient.fetch(qr1)
    assert len(res) == len(qr1)

    res2 = LCClient.fetch(qr1[0])
    assert len(res2) == 1


@pytest.mark.remote_data
@pytest.mark.parametrize(
    'query',
    [(a.Time('2012/10/4', '2012/10/5') & a.Instrument.eve & a.Level.zero)])
def test_fido(LCClient, query):
    qr = Fido.search(query)
    client = qr[0].client
    assert isinstance(qr, UnifiedResponse)
    assert isinstance(client, eve.EVEClient)
    response = Fido.fetch(qr)
    assert len(response) == qr._numfile


@pytest.mark.remote_data
@pytest.mark.parametrize(
    'time',
    [(a.Time('2012/10/4', '2012/10/6'))])
def test_levels(time):
    """
    Test the correct handling of level
    Level 0 comes from EVEClient, other levels from EVE.
    """
    eve_a = a.Instrument.eve
    qr = Fido.search(time, eve_a, a.Level.one)
    clients = {type(a.client) for a in qr}
    assert clients == {VSOClient}

    qr = Fido.search(time, eve_a, a.Level.zero)
    clients = {type(a.client) for a in qr}
    assert clients == {eve.EVEClient}


def test_attr_reg():
    assert a.Instrument.eve == a.Instrument('EVE')
    assert a.Level.zero == a.Level('0')


def test_client_repr(LCClient):
    """
    Repr check
    """
    output = str(LCClient)
    assert output[:50] == 'sunpy.net.dataretriever.sources.eve.EVEClient\n\nPro'


def mock_query_object(LCClient):
    """
    Creating a Query Response object and prefilling it with some information
    """
    # Creating a Query Response Object
    start = '2016/1/1'
    end = '2016/1/1 23:59:59'
    obj = {
        'Start Time': parse_time(start),
        'End Time': parse_time(end),
        'Instrument': 'EVE',
        'Physobs': 'irradiance',
        'Source': 'SDO',
        'Provider': 'LASP',
        'Level': '0',
        'url': ('http://lasp.colorado.edu/eve/data_access/evewebdata/'
                'quicklook/L0CS/SpWx/2016/20160101_EVE_L0CS_DIODES_1m.txt')
    }
    results = QueryResponse([obj], client=LCClient)
    return results


def test_show(LCClient):
    mock_qr = mock_query_object(LCClient)
    qrshow0 = mock_qr.show()
    qrshow1 = mock_qr.show('Start Time', 'Instrument')
    allcols = {'Start Time', 'End Time', 'Instrument', 'Physobs', 'Source',
               'Provider', 'Level', 'url'}
    assert not allcols.difference(qrshow0.colnames)
    assert qrshow1.colnames == ['Start Time', 'Instrument']
    assert qrshow0['Instrument'][0] == 'EVE'
