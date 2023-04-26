import pytest
from hypothesis import given

import astropy.units as u
from astropy.time import TimeDelta

import sunpy.net.dataretriever.sources.lyra as lyra
from sunpy.net import Fido
from sunpy.net import attrs as a
from sunpy.net._attrs import Instrument, Time
from sunpy.net.dataretriever.client import QueryResponse
from sunpy.net.fido_factory import UnifiedResponse
from sunpy.net.tests.strategies import range_time
from sunpy.time import parse_time


@pytest.fixture
def LCClient():
    return lyra.LYRAClient()


@pytest.mark.remote_data
@pytest.mark.parametrize(("timerange", "url_start", "url_end"), [
    (Time('2012/1/7', '2012/1/7'),
     'http://proba2.oma.be/lyra/data/bsd/2012/01/07/lyra_20120107-000000_lev2_std.fits',
     'http://proba2.oma.be/lyra/data/bsd/2012/01/07/lyra_20120107-000000_lev2_std.fits'
     ),
    (Time('2012/12/1', '2012/12/2'),
     'http://proba2.oma.be/lyra/data/bsd/2012/12/01/lyra_20121201-000000_lev2_std.fits',
     'http://proba2.oma.be/lyra/data/bsd/2012/12/02/lyra_20121202-000000_lev2_std.fits'
     ),
    (Time('2012/4/7', '2012/4/14'),
     'http://proba2.oma.be/lyra/data/bsd/2012/04/07/lyra_20120407-000000_lev2_std.fits',
     'http://proba2.oma.be/lyra/data/bsd/2012/04/14/lyra_20120414-000000_lev2_std.fits'
     )
])
def test_get_url_for_time_range(LCClient, timerange, url_start, url_end):
    qresponse = LCClient.search(timerange, a.Level.two)
    urls = [i['url'] for i in qresponse]
    assert isinstance(urls, list)
    assert urls[0] == url_start
    assert urls[-1] == url_end


@given(range_time('2010-01-06'))
def test_can_handle_query(time):
    LCClient = lyra.LYRAClient()
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
    assert qr1[0]['Start Time'] == time.start
    almost_day = TimeDelta(1 * u.day - 1 * u.millisecond)
    assert qr1[-1]['End Time'] == time.end + almost_day


@pytest.mark.remote_data
@pytest.mark.parametrize(("time", "instrument"), [
    (Time('2013/8/27', '2013/8/27'), Instrument('lyra'))])
def test_get(LCClient, time, instrument):
    # Cut it down to one file with the level
    qr1 = LCClient.search(time, instrument, a.Level.two)
    download_list = LCClient.fetch(qr1)
    assert len(download_list) == len(qr1)


@pytest.mark.remote_data
@pytest.mark.parametrize(
    ("time", "instrument"),
    [(a.Time('2012/10/4', '2012/10/4'), a.Instrument.lyra)])
def test_fido(time, instrument):
    # Cut it down to one file with the level
    qr = Fido.search(time, instrument, a.Level.two)
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
        'Instrument': 'LYRA',
        'Physobs': 'irradiance',
        'Source': 'PROBA2',
        'Provider': 'ESA',
        'Level': '2',
        'url': ('http://proba2.oma.be/lyra/data/bsd/2016/01/01/'
                'lyra_20160101-000000_lev2_std.fits')
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
    assert qrshow0['Instrument'][0] == 'LYRA'
