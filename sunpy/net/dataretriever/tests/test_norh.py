from unittest import mock
from datetime import datetime, timedelta

import pytest
from hypothesis import given, settings, HealthCheck

import astropy.units as u
from astropy.time import TimeDelta
from astropy.time import Time
import sunpy.net.dataretriever.sources.norh as norh
from sunpy.net.download import Results
from sunpy.net.tests.strategies import time_attr, range_time
from sunpy.net.dataretriever.client import QueryResponse
from sunpy.net.fido_factory import UnifiedResponse
from sunpy.net import Fido
from sunpy.net import attrs as a
from sunpy.time import parse_time
from sunpy.time.timerange import TimeRange

NORHClient = norh.NoRHClient()
BASEURL = 'ftp://anonymous:data@sunpy.org@solar-pub.nao.ac.jp/pub/nsro/norh/data/tcx/%Y/%m/{freq}%y%m%d'

def create_url(start, end, wavelength):
    """
    This function creates a url based on the NoRHClient data,
    instead of making an online request.
    """
    if wavelength == 34 * u.GHz:
            freq = 'tcz'
    elif wavelength == 17 * u.GHz:
            freq = 'tca'
    else:
        raise ValueError('Wavelength should be 17Ghz or 34Ghz, found {}'.format(wavelength))
    start = datetime.strptime(start, '%Y/%m/%d')
    end = datetime.strptime(end, '%Y/%m/%d')

    lst = list()
    for n in range(int((end - start).days)+1):
        lst.append((start+timedelta(n)).strftime(BASEURL.format(freq=freq)))
    return lst


def mock_query_object(start_date, end_date, wavelength):
    """
    Creation of a QueryResponse object, and prefill some
    downloaded data from norh.NoRHClient.fetch(Time('20 ..))
    """
    # Create a mock QueryResponse object
    map_ = {
        'TimeRange': TimeRange(parse_time(start_date), parse_time(end_date)),
        'Time_start': parse_time(start_date),
        'Time_end':  parse_time(end_date),
        'source': 'NAOJ',
        'instrument': 'NORH',
        'physobs': '',
        'provider': 'NRO',
        'wavelength': wavelength,
    }

    # Create a range of Time
    st_datetime = datetime.strptime(start_date, '%Y/%m/%d')
    ed_datetime = datetime.strptime(end_date, '%Y/%m/%d') + timedelta(days=1)
    time_range = TimeRange(st_datetime, ed_datetime).split(int((ed_datetime-st_datetime).days))

    with mock.patch('sunpy.net.dataretriever.sources.norh.NoRHClient._get_url_for_timerange',
                    return_value=create_url(start_date, end_date, wavelength)):
        resp = QueryResponse.create(map_,
                                    NORHClient._get_url_for_timerange(
                                        map_['TimeRange'], wavelength=wavelength),
                                    time=time_range)
        # Attach the client with the QueryResponse
        resp.client = NORHClient
        return resp


@pytest.mark.remote_data
def test_fetch_working():
    """
    Tests if the mock fetch contains the correct data.
    """

    qr1 = NORHClient.search(a.Time('2012/10/4', '2012/10/6'),
                            a.Instrument('norh'), a.Wavelength(17*u.GHz))

    # Create a mock query object
    mock_qr = mock_query_object('2012/10/4', '2012/10/6', wavelength=17*u.GHz)

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
    res = NORHClient.fetch(qr1)
    download_list = res.wait(progress=False)
    assert len(download_list) == len(qr1)

    download_list.sort()
    # Assert each of the files have the same file names
    assert download_list[0].split('/')[-1] == 'tca121004'
    assert download_list[1].split('/')[-1] == 'tca121005'
    assert download_list[2].split('/')[-1] == 'tca121006'

@pytest.fixture
def create_mock_url(mocker,sdate,edate,wave):
    mocker.patch('sunpy.net.dataretriever.sources.norh.NoRHClient._get_url_for_timerange',
                return_value=create_url(sdate,edate,wave))

@pytest.fixture
def create_mock_fetch(mocker, sdate, edate, wave):
    mocker.patch('sunpy.net.fido_factory.Fido.fetch', side_effect=(
           UnifiedResponse(mock_query_object(sdate, edate, wave))))

@pytest.mark.usefixtures('create_mock_url')
@pytest.mark.parametrize(
    "sdate, edate, wave",
    [('2012/3/7', '2012/3/14',17*u.GHz),
    ('2012/3/7', '2012/3/14', 34*u.GHz)]
)
def test_get_url_for_time_range(sdate,edate,wave):
    urls = norh.NoRHClient()._get_url_for_timerange(TimeRange(parse_time('2012/3/7'),
                                                              parse_time('2012/3/14')),
                                                    wavelength=wave)
    assert isinstance(urls, list)
    if wave.value == 17.0:
        assert urls[0] == 'ftp://anonymous:data@sunpy.org@solar-pub.nao.ac.jp/pub/nsro/norh/data/tcx/2012/03/tca120307'
        assert urls[-1] == 'ftp://anonymous:data@sunpy.org@solar-pub.nao.ac.jp/pub/nsro/norh/data/tcx/2012/03/tca120314'
    else:
        assert urls[0] == 'ftp://anonymous:data@sunpy.org@solar-pub.nao.ac.jp/pub/nsro/norh/data/tcx/2012/03/tcz120307'
        assert urls[-1] == 'ftp://anonymous:data@sunpy.org@solar-pub.nao.ac.jp/pub/nsro/norh/data/tcx/2012/03/tcz120314'


@given(time_attr())
def test_can_handle_query(time):
    ans1 = norh.NoRHClient._can_handle_query(time, a.Instrument('norh'))
    assert ans1 is True
    ans1 = norh.NoRHClient._can_handle_query(time, a.Instrument('norh'),
                                             a.Wavelength(10*u.GHz))
    assert ans1 is True
    ans2 = norh.NoRHClient._can_handle_query(time)
    assert ans2 is False

@pytest.mark.remote_data
@pytest.mark.parametrize("wave", [a.Wavelength(17*u.GHz), a.Wavelength(34*u.GHz)])
@given(time=range_time(Time('1992-6-1')))
@settings(max_examples=2, deadline=50000)
def test_query(time, wave):
    qr1 = norh.NoRHClient().search(time, a.Instrument('norh'), wave)
    assert isinstance(qr1, QueryResponse)
    # Not all hypothesis queries are going to produce results, and
    if qr1:
        # There are no observations everyday
        #  so the results found have to be equal or later than the queried time
        #  (looking at the date because it may search for miliseconds, but only date is available)
        assert qr1.time_range().start.strftime('%Y-%m-%d') >= time.start.strftime('%Y-%m-%d')
        #  and the end time equal or smaller.
        # hypothesis can give same start-end, but the query will give you from start to end (so +1)
        assert qr1.time_range().end <= time.end + TimeDelta(1*u.day)


# Don't use time_attr here for speed.
def test_query_no_wave():
    c = norh.NoRHClient()
    with pytest.raises(ValueError):
        c.search(a.Time("2016/10/1", "2016/10/2"), a.Instrument('norh'))


def test_wavelength_range():
    with pytest.raises(ValueError):
        norh.NoRHClient().search(
            a.Time("2016/10/1", "2016/10/2"), a.Instrument('norh'),
            a.Wavelength(17 * u.GHz, 34 * u.GHz))


def test_query_wrong_wave():
    c = norh.NoRHClient()
    with pytest.raises(ValueError):
        c.search(a.Time("2016/10/1", "2016/10/2"), a.Instrument('norh'), a.Wavelength(50*u.GHz))


@mock.patch('sunpy.net.dataretriever.sources.norh.NoRHClient._get_url_for_timerange',
            side_effect=create_url('2013/10/5', '2013/10/7', wavelength=34*u.GHz))
@mock.patch('sunpy.net.dataretriever.sources.norh.NoRHClient.search',
            side_effect=(UnifiedResponse(mock_query_object('2013/10/5', '2013/10/7',
                         wavelength=34*u.GHz))))
@mock.patch('sunpy.net.download.Results.wait', return_value=['some/path/extension/tcz131005',
            'some/path/extension/tcz131007', 'some/path/extension/tcz131006'])
def test_get(mock_result, mock_fetch, mock_timerange):
    LCClient = norh.NoRHClient()
    qr1 = LCClient.search(a.Time('2013/10/5', '2013/10/7'), a.Instrument('norh'),
                          a.Wavelength(34*u.GHz))
    res = LCClient.fetch(qr1)
    assert isinstance(res, Results)
    download_list = res.wait(progress=False)
    assert len(download_list) == len(qr1)


@pytest.mark.usefixtures('create_mock_fetch')
@pytest.mark.usefixtures('create_mock_url')
@pytest.mark.parametrize(
    "sdate, edate, wave",
    [('2012/10/4', '2012/10/6',17*u.GHz),
    ('2012/10/4', '2012/10/6', 34*u.GHz)]
)
def test_fido_34(sdate,edate, wave):
    qr = Fido.search(a.Time('2012/10/4', '2012/10/6'), a.Instrument('norh'), a.Wavelength(wave))
    assert isinstance(qr, UnifiedResponse)
    response = Fido.fetch(qr)
    assert len(response) == qr._numfile
