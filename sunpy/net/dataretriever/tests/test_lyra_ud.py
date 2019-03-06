import pytest
from unittest import mock
from datetime import datetime, timedelta

import astropy.units as u
from sunpy.time.timerange import TimeRange
from sunpy.time import parse_time
from sunpy.net.vso.attrs import Time, Instrument
from sunpy.net.dataretriever.client import QueryResponse
import sunpy.net.dataretriever.sources.lyra as lyra
from sunpy.net.fido_factory import UnifiedResponse
from sunpy.net import Fido
from sunpy.net import attrs as a

from hypothesis import given, settings
from sunpy.net.tests.strategies import time_attr

LCClient = lyra.LYRAClient()


def mock_query_object(start_date, end_date):
    """
    Creation of a QueryResponse object, and prefill some
    downloaded data from lyra.LYRAClient().fetch(Time('20 ..)
    """
    # Create a mock QueryResponse object
    map_ = {
        'TimeRange': TimeRange(parse_time(start_date), parse_time(end_date)),
        'Time_start': parse_time(start_date),
        'Time_end':  parse_time(end_date),
        'source': 'Proba2',
        'instrument': 'lyra',
        'physobs': 'irradiance',
        'provider': 'esa'
    }

    resp = QueryResponse.create(map_,
    LCClient._get_url_for_timerange(TimeRange(parse_time(start_date), parse_time(end_date))))
    # Attach the client with the QueryResponse
    resp.client = LCClient
    return resp


@pytest.mark.remote_data
def test_fetch_working():
    """
    Tests if the mock fetch contains the correct data.
    """

    qr1 = LCClient.search(a.Time('2012/10/4', '2012/10/6'),
                          a.Instrument('lyra'))

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
    assert qr1.time_range() == TimeRange('2012/10/4', '2012/10/6')

    # Assert the fetch object, and whether it returns the correct set of files
    res = LCClient.fetch(qr1)
    download_list = res.wait(progress=False)
    assert len(download_list) == len(qr1)


def daterange(start_date, end_date):
    """
    Helper function to create an iterable of datetime
    containing each days as its entries from start to end
    date
    """

    start_date = datetime.strptime(start_date, '%Y/%m/%d')
    end_date = datetime.strptime(end_date, '%Y/%m/%d')
    range_of_date = int((end_date - start_date).days)

    if range_of_date is not 0:
        for n in range(range_of_date+1):
            yield start_date + timedelta(n)
    else:
        yield start_date


# Create mocker fixtures for using parametrized tests
@pytest.fixture
def create_mock_fetch(mocker, sdate, edate, instrument):
    mocker.patch('sunpy.net.fido_factory.Fido.fetch',
        side_effect=(UnifiedResponse(mock_query_object(sdate, edate))))


@pytest.fixture
def create_mock_search(mocker, sdate, edate, instrument):
    mocker.patch('sunpy.net.dataretriever.sources.lyra.LYRAClient.search',
        return_value=mock_query_object(sdate, edate))


@pytest.fixture
def create_mock_result(mocker, sdate, edate, instrument):
    mocker.patch('sunpy.net.download.Results.wait',
    return_value=[LCClient._get_url_for_date(parse_time(date)) for date in daterange(sdate, edate)])


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
def test_get_url_for_time_range(timerange, url_start, url_end):
    urls = LCClient._get_url_for_timerange(timerange)
    assert isinstance(urls, list)
    assert urls[0] == url_start
    assert urls[-1] == url_end


def test_get_url_for_date():
    url = LCClient._get_url_for_date(parse_time((2013, 2, 13)))
    assert url == 'http://proba2.oma.be/lyra/data/bsd/2013/02/13/lyra_20130213-000000_lev2_std.fits'


@given(time_attr())
def test_can_handle_query(time):
    ans1 = lyra.LYRAClient._can_handle_query(
        time, Instrument('lyra'))
    assert ans1 is True
    ans2 = lyra.LYRAClient._can_handle_query(time)
    assert ans2 is False


@settings(deadline=50000)
@given(time_attr())
def test_query(time):
    qr1 = LCClient.search(time, Instrument('lyra'))
    assert isinstance(qr1, QueryResponse)
    assert qr1.time_range().start == time.start
    assert qr1.time_range().end == time.end


@pytest.mark.usefixtures('create_mock_search')
@pytest.mark.usefixtures('create_mock_result')
@pytest.mark.parametrize("sdate, edate ,instrument", [
    ('2013/8/27', '2013/8/27', Instrument('lyra')),
    ('2013/2/4', '2013/2/6', Instrument('lyra')),
])
def test_get(sdate, edate, instrument):
    qr1 = LCClient.search(TimeRange(sdate, edate), instrument)
    res = LCClient.fetch(qr1)
    download_list = res.wait(progress=False)
    assert len(download_list) == len(qr1)


@pytest.mark.usefixtures('create_mock_fetch')
@pytest.mark.parametrize(
    "sdate, edate, instrument",
    [('2012/10/4', '2012/10/6', a.Instrument('lyra')),
     ('2013/10/5', '2013/10/7', a.Instrument('lyra'))]
     )
def test_fido(sdate, edate, instrument):
    qr = Fido.search(a.Time(sdate, edate), a.Instrument('lyra'))
    assert isinstance(qr, UnifiedResponse)
    response = Fido.fetch(qr)
    assert len(response) == qr._numfile
