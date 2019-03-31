import pytest
from unittest import mock
from datetime import datetime, timedelta

from sunpy.time import parse_time
from sunpy.time.timerange import TimeRange
from sunpy.net.vso.vso import VSOClient
from sunpy.net.vso.attrs import Time, Instrument, Source, Level
from sunpy.net.dataretriever.client import QueryResponse
import sunpy.net.dataretriever.sources.eve as eve
from sunpy.net.fido_factory import UnifiedResponse
from sunpy.net import Fido
from sunpy.net import attr
from sunpy.net import attrs as a


LCClient = eve.EVEClient()
BASEURL = eve.BASEURL


def create_url(start, end):
    """
    This function creates a url based on the EVEClient data,
    instead of making an online request.
    """

    start = datetime.strptime(start, '%Y/%m/%d')
    end = datetime.strptime(end, '%Y/%m/%d')

    lst = list()
    for n in range(int((end - start).days)+1):
        lst.append((start+timedelta(n)).strftime(BASEURL))
    return lst


def mock_query_object(start_date, end_date, level=0):
    """
    Creation of a QueryResponse object, and prefill some
    downloaded data from eve.EVEClient().fetch(Time('20 ..)
    """
    # Create a mock QueryResponse object
    map_ = {
        'TimeRange': TimeRange(start_date, end_date),
        'Time_start': start_date,
        'Time_end':  end_date,
        'source': 'SDO',
        'instrument': 'eve',
        'physobs': 'irradiance',
        'provider': 'LASP'
    }
    st_datetime = datetime.strptime(start_date, '%Y/%m/%d')
    ed_datetime = datetime.strptime(end_date, '%Y/%m/%d') + timedelta(days=1)
    time_range = TimeRange(st_datetime, ed_datetime).split(int((ed_datetime-st_datetime).days))

    with mock.patch('sunpy.net.dataretriever.sources.eve.EVEClient._get_url_for_timerange',
                    return_value=(create_url(start_date, end_date))):
        resp = QueryResponse.create(map_,
        LCClient._get_url_for_timerange(TimeRange(start_date, end_date)),
        time=time_range)
    # Attach the client with the QueryResponse
    # Here level 0 corresponds to the default client, i.e. LCClient
    if level is 0:
        resp.client = LCClient
    else:
        resp.client = VSOClient()
    return resp


@pytest.mark.remote_data
def test_fetch_working():
    """
    Tests if the mock fetch contains the correct data.
    """

    qr1 = LCClient.search(a.Time('2012/10/4', '2012/10/6'),
                            a.Instrument('eve'))

    # Create a mock query object
    mock_qr = mock_query_object('2012/10/4', '2012/10/6')

    mock_qr = mock_qr[0]
    qr = qr1[0]
    # Assert the values
    assert mock_qr.source == qr.source
    assert mock_qr.instrument == qr.instrument
    assert mock_qr.physobs == qr.physobs
    assert mock_qr.provider == qr.provider
    assert mock_qr.url == qr.url
    assert mock_qr.time == qr.time

    # Assert if the time range is same
    assert qr1.time_range() == TimeRange('2012/10/4', '2012/10/7')

    # Assert the fetch object, and whether it returns the correct set of files
    res = LCClient.fetch(qr1)
    download_list = res.wait(progress=False)
    assert len(download_list) == len(qr1)


@pytest.fixture
def create_mock_url(mocker, sdate, edate, url_start, url_end):
    mocker.patch('sunpy.net.dataretriever.sources.eve.EVEClient._get_url_for_timerange',
                    return_value=(create_url(sdate, edate)))


@pytest.fixture
def create_mock_search(mocker, sdate, edate):
    mocker.patch('sunpy.net.dataretriever.sources.eve.EVEClient.search',
            return_value=mock_query_object(sdate, edate))


@pytest.fixture
def create_mock_search_using_side_effect(mocker, sdate, edate):
    mocker.patch('sunpy.net.Fido.search',
            side_effect=sid_effect)


@pytest.mark.usefixtures('create_mock_url')
@pytest.mark.parametrize("sdate, edate, url_start, url_end", [
    ('2012/4/21', '2012/4/21',
     'http://lasp.colorado.edu/eve/data_access/evewebdata/quicklook/L0CS/SpWx/2012/20120421_EVE_L0CS_DIODES_1m.txt',
     'http://lasp.colorado.edu/eve/data_access/evewebdata/quicklook/L0CS/SpWx/2012/20120421_EVE_L0CS_DIODES_1m.txt'
     ),
    ('2012/5/5', '2012/5/6',
     'http://lasp.colorado.edu/eve/data_access/evewebdata/quicklook/L0CS/SpWx/2012/20120505_EVE_L0CS_DIODES_1m.txt',
     'http://lasp.colorado.edu/eve/data_access/evewebdata/quicklook/L0CS/SpWx/2012/20120506_EVE_L0CS_DIODES_1m.txt',
     ),
    ('2012/7/7', '2012/7/14',
     'http://lasp.colorado.edu/eve/data_access/evewebdata/quicklook/L0CS/SpWx/2012/20120707_EVE_L0CS_DIODES_1m.txt',
     'http://lasp.colorado.edu/eve/data_access/evewebdata/quicklook/L0CS/SpWx/2012/20120714_EVE_L0CS_DIODES_1m.txt',
     )
])
def test_get_url_for_time_range(sdate, edate, url_start, url_end):
    urls = LCClient._get_url_for_timerange(TimeRange(sdate, edate))

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


@mock.patch('sunpy.net.dataretriever.sources.eve.EVEClient.search',
            return_value=mock_query_object('2012/8/9', '2012/8/10'))
def test_query(mock_search):
    qr1 = LCClient.search(Time('2012/8/9', '2012/8/10'), Instrument('eve'))
    assert isinstance(qr1, QueryResponse)
    assert len(qr1) == 2
    assert qr1.time_range().start == parse_time('2012/08/09')
    assert qr1.time_range().end == parse_time('2012/08/11')  # includes end.



@mock.patch('sunpy.net.dataretriever.sources.eve.EVEClient.search',
            return_value=mock_query_object('2012/11/27', '2012/11/27'))
@mock.patch('sunpy.net.download.Results.wait',
            return_value=['some/path/extension/20121127_EVE_L0CS_DIODES_1m.1.txt'])
def test_get(time, instrument):
    qr1 = LCClient.search(time, instrument)
    res = LCClient.fetch(qr1)
    download_list = res.wait(progress=False)
    assert len(download_list) == len(qr1)


@mock.patch('sunpy.net.dataretriever.sources.eve.EVEClient._get_url_for_timerange',
            return_value=(create_url('2012/10/4', '2012/10/6')))
@mock.patch('sunpy.net.fido_factory.Fido.search', return_value=(
           UnifiedResponse(mock_query_object('2012/10/4', '2012/10/6'))))
@mock.patch('sunpy.net.download.Results.wait',
            return_value=['some/path/extension/20121006_EVE_L0CS_DIODES_1m.1.txt',
             'some/path/extension/20121005_EVE_L0CS_DIODES_1m.1.txt',
             'some/path/extension/20121004_EVE_L0CS_DIODES_1m.1.txt']
)
def test_fido(mock_result, mock_search, mock_url):
    qr = Fido.search((a.Time('2012/10/4', '2012/10/6') & a.Instrument('eve') & a.Level(0)))
    client = qr.get_response(0).client
    assert isinstance(qr, UnifiedResponse)
    assert isinstance(client, eve.EVEClient)
    response = Fido.fetch(qr)
    assert len(response) == qr._numfile


def sid_effect(time, *args):
    """
    Helper function to return a list of UnifiedResponse object based on the
    `sunpy.net.attrs.Level`.
    """

    level = args[-1]
    LevelVal = None

    if isinstance(level, attr.AttrOr):
        LevelVal = [level.attrs[0].value, level.attrs[1].value]
        qr = [mock_query_object(time[0], time[1], LevelVal[0]),
                 mock_query_object(time[0], time[1], LevelVal[1])]

    else:
        LevelVal = level.value
        qr = mock_query_object(time[0], time[1], LevelVal)

    return UnifiedResponse(qr)


@pytest.mark.usefixtures('create_mock_url')
@pytest.mark.usefixtures('create_mock_search_using_side_effect')
@pytest.mark.parametrize(
    'sdate, edate, url_start, url_end',
    [('2012/10/4', '2012/10/6', '', ''), ('2012/11/27', '2012/11/27', '', '')])
def test_levels(sdate, edate, url_start, url_end):
    """
    Test the correct handling of level 0 / 1.
    The default should be level 1 from VSO, level 0 comes from EVEClient.
    """
    time = sdate, edate
    eve_a = a.Instrument('EVE')
    qr = Fido.search(time, eve_a, a.Level(1))

    client = qr.get_response(0).client
    assert isinstance(client, VSOClient)

    qr = Fido.search(time, eve_a, a.Level(0))
    client = qr.get_response(0).client
    assert isinstance(client, eve.EVEClient)

    qr = Fido.search(time, eve_a, a.Level(0) | a.Level(1))
    clients = {type(a.client) for a in qr.responses}
    assert clients.symmetric_difference({VSOClient, eve.EVEClient}) == set()
