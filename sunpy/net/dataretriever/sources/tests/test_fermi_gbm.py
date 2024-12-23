import pytest
from hypothesis import given

import astropy.units as u
from astropy.time import TimeDelta

import sunpy.net.dataretriever.sources.fermi_gbm as fermi_gbm
from sunpy.net import Fido
from sunpy.net import attrs as a
from sunpy.net.dataretriever.client import QueryResponse
from sunpy.net.fido_factory import UnifiedResponse
from sunpy.net.tests.strategies import time_attr
from sunpy.time import parse_time


@pytest.fixture
def LCClient():
    return fermi_gbm.GBMClient()


@pytest.mark.remote_data
@pytest.mark.parametrize(("timerange", "url_start", "url_end"),
                         [(a.Time('2011/06/07', '2011/06/09'),
                           'https://heasarc.gsfc.nasa.gov/FTP/fermi/data/gbm/daily/2011/06/07/'
                           'current/glg_cspec_n5_110607_v00.pha',
                           'https://heasarc.gsfc.nasa.gov/FTP/fermi/data/gbm/daily/2011/06/09/'
                           'current/glg_cspec_n5_110609_v00.pha')])
def test_get_url_for_time_range(LCClient, timerange, url_start, url_end):
    qresponse = LCClient.search(timerange, a.Detector.n5, a.Resolution.cspec)
    urls = [i['url'] for i in qresponse]
    assert isinstance(urls, list)
    assert urls[0] == url_start
    assert urls[-1] == url_end


@given(time_attr())
def test_can_handle_query(time):
    LCClient = fermi_gbm.GBMClient()
    ans1 = LCClient._can_handle_query(time, a.Instrument.gbm)
    assert ans1 is True
    ans2 = LCClient._can_handle_query(time, a.Instrument.gbm,
                                      a.Detector.n5)
    assert ans2 is True
    ans3 = LCClient._can_handle_query(time, a.Instrument.gbm,
                                      a.Detector.n5, a.Resolution.ctime)
    assert ans3 is True
    ans4 = LCClient._can_handle_query(time)
    assert ans4 is False


@pytest.mark.remote_data
@pytest.mark.parametrize(("time", "instrument"), [
    (a.Time('2012/8/9', '2012/8/9'), a.Instrument.gbm),
])
def test_query(LCClient, time, instrument):
    qr1 = LCClient.search(time, instrument, a.Detector.n5, a.Resolution.ctime)
    assert isinstance(qr1, QueryResponse)
    assert len(qr1) == 1
    almost_day = TimeDelta(1 * u.day - 1 * u.millisecond)
    assert qr1.time_range().start == time.start.to_datetime()
    assert qr1.time_range().end == (time.end + almost_day).to_datetime()


@pytest.mark.remote_data
@pytest.mark.parametrize(("time", "instrument"), [
    (a.Time('2012/11/27', '2012/11/27'), a.Instrument.gbm),
])
def test_get(LCClient, time, instrument):
    qr1 = LCClient.search(time, instrument, a.Detector.n5, a.Resolution.ctime)
    download_list = LCClient.fetch(qr1)
    assert len(download_list) == len(qr1)


@pytest.mark.remote_data
@pytest.mark.parametrize(
    'query',
    [(a.Time('2012/10/4', '2012/10/5') & a.Instrument.gbm & a.Detector.n5)])
def test_fido(LCClient, query):
    qr = Fido.search(query)
    client = qr[0].client
    assert isinstance(qr, UnifiedResponse)
    assert isinstance(client, type(LCClient))
    response = Fido.fetch(qr)
    assert len(response) == qr._numfile


def test_attr_reg():
    assert a.Instrument.gbm == a.Instrument('GBM')


def test_client_repr(LCClient):
    """
    Repr check
    """
    output = str(LCClient)
    assert output[:51] == 'sunpy.net.dataretriever.sources.fermi_gbm.GBMClient'


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
        'Instrument': 'GBM',
        'Physobs': 'flux',
        'Source': 'FERMI',
        'Provider': 'NASA',
        'Resolution': 'cspec',
        'Detector': 'n5',
        'url': ('https://heasarc.gsfc.nasa.gov/FTP/fermi/data/gbm/daily/'
                '2016/01/01/current/glg_cspec_n5_160101_v00.pha')
    }
    results = QueryResponse([obj], client=LCClient)
    return results


def test_show(LCClient):
    mock_qr = mock_query_object(LCClient)
    qrshow0 = mock_qr.show()
    qrshow1 = mock_qr.show('Start Time', 'Instrument')
    allcols = {'Start Time', 'End Time', 'Instrument', 'Physobs', 'Source',
               'Provider', 'Resolution', 'Detector', 'url'}
    assert not allcols.difference(qrshow0.colnames)
    assert qrshow1.colnames == ['Start Time', 'Instrument']
    assert qrshow0['Instrument'][0] == 'GBM'
