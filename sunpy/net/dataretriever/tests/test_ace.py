import pytest

import astropy.units as u
from astropy.time import TimeDelta

from sunpy.time.timerange import TimeRange
from sunpy.time import parse_time
from sunpy.net.vso.attrs import Time, Instrument
from sunpy.net.dataretriever.client import QueryResponse
from sunpy.net.fido_factory import UnifiedResponse
from sunpy.net import Fido, attrs as a

import sunpy.net.dataretriever.sources.ace as ace


@pytest.fixture
def SWEPAMClient():
    return ace.SWEPAMClient()


@pytest.fixture
def EPAMClient():
    return ace.EPAMClient()


@pytest.fixture
def MAGClient():
    return ace.MAGClient()


@pytest.fixture
def SISClient():
    return ace.SISClient()


@pytest.fixture
def client(request):
    return request.getfixturevalue(request.param)


@pytest.mark.remote_data
@pytest.mark.parametrize("timerange, url_start, url_end",
                         [(TimeRange('2015/12/27', '2015/12/30'),
                           'ftp://ftp.swpc.noaa.gov/pub/lists/ace/20151227_ace_swepam_1m.txt',
                           'ftp://ftp.swpc.noaa.gov/pub/lists/ace/20151230_ace_swepam_1m.txt')])
def test_get_url_for_timerange(SWEPAMClient, timerange, url_start, url_end):
    urls = SWEPAMClient._get_url_for_timerange(timerange)
    assert isinstance(urls, list)
    assert urls[0] == url_start
    assert urls[-1] == url_end


TRANGE = Time('2015/12/30', '2015/12/31')


@pytest.mark.parametrize("client, time, instrument, expected",
                         [('SWEPAMClient', TRANGE, Instrument('swepam'), True),
                          ('EPAMClient', TRANGE, Instrument('epam'), True),
                          ('MAGClient', TRANGE, Instrument('mag'), True),
                          ('SISClient', TRANGE, Instrument('sis'), True),
                          ('SWEPAMClient', TRANGE, Instrument('swap'), False),
                          ('EPAMClient', TRANGE, None, False)], indirect=["client"]
                         )
def test_can_handle_query(client, time, instrument, expected):
    assert client._can_handle_query(time, instrument) is expected


@pytest.mark.remote_data
def test_query(SWEPAMClient):
    qr = SWEPAMClient.search(Time('2015/12/27', '2015/12/30'), a.Instrument('swepam'))
    assert isinstance(qr, QueryResponse)
    assert len(qr) == 4
    assert qr.time_range().start == parse_time('2015/12/27')
    almost_day = TimeDelta(1 * u.day - 1 * u.millisecond)
    assert qr.time_range().end == parse_time('2015/12/30') + almost_day


@pytest.mark.remote_data
@pytest.mark.parametrize("time, instrument",
                         [(Time('2015/12/27', '2015/12/30'), a.Instrument('swepam'))])
def test_get_swepam(SWEPAMClient, time, instrument):
    qr = SWEPAMClient.search(time, instrument)
    download_list = SWEPAMClient.fetch(qr)
    assert len(download_list) == len(qr)


@pytest.mark.remote_data
def test_fido_query_swepam(SWEPAMClient):
    qr = Fido.search(a.Time('2015/12/27', '2015/12/30'), a.Instrument('swepam'))
    assert isinstance(qr, UnifiedResponse)
    response = Fido.fetch(qr)
    assert len(response) == qr._numfile


@pytest.mark.remote_data
@pytest.mark.parametrize("time, instrument",
                         [(Time('2015/12/27', '2015/12/30'), a.Instrument('epam'))])
def test_get_epam(EPAMClient, time, instrument):
    qr = EPAMClient.search(time, instrument)
    download_list = EPAMClient.fetch(qr)
    assert len(download_list) == len(qr)


@pytest.mark.remote_data
def test_fido_query_epam(EPAMClient):
    qr = Fido.search(a.Time('2015/12/27', '2015/12/30'), a.Instrument('epam'))
    assert isinstance(qr, UnifiedResponse)
    response = Fido.fetch(qr)
    assert len(response) == qr._numfile


@pytest.mark.remote_data
@pytest.mark.parametrize("time, instrument",
                         [(Time('2015/12/27', '2015/12/30'), a.Instrument('mag'))])
def test_get_mag(MAGClient, time, instrument):
    qr = MAGClient.search(time, instrument)
    download_list = MAGClient.fetch(qr)
    assert len(download_list) == len(qr)


@pytest.mark.remote_data
def test_fido_query_mag(MAGClient):
    qr = Fido.search(a.Time('2015/12/27', '2015/12/30'), a.Instrument('mag'))
    assert isinstance(qr, UnifiedResponse)
    response = Fido.fetch(qr)
    assert len(response) == qr._numfile


@pytest.mark.remote_data
@pytest.mark.parametrize("time, instrument",
                         [(Time('2015/12/27', '2015/12/30'), a.Instrument('sis'))])
def test_get_sis(SISClient, time, instrument):
    qr = SISClient.search(time, instrument)
    download_list = SISClient.fetch(qr)
    assert len(download_list) == len(qr)


@pytest.mark.remote_data
def test_fido_query_sis(SISClient):
    qr = Fido.search(a.Time('2015/12/27', '2015/12/30'), a.Instrument('sis'))
    assert isinstance(qr, UnifiedResponse)
    response = Fido.fetch(qr)
    assert len(response) == qr._numfile
