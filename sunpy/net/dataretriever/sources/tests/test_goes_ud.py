import pytest
from hypothesis import given

import astropy.units as u
from astropy.time import TimeDelta

import sunpy.net.dataretriever.sources.goes as goes
from sunpy.net import Fido
from sunpy.net import attrs as a
from sunpy.net._attrs import Instrument, Time
from sunpy.net.dataretriever.client import QueryResponse
from sunpy.net.fido_factory import UnifiedResponse
from sunpy.net.tests.strategies import goes_time
from sunpy.time import is_time_equal, parse_time
from sunpy.time.timerange import TimeRange


@pytest.fixture
def LCClient():
    return goes.XRSClient()


@pytest.mark.remote_data
@pytest.mark.parametrize(
    "timerange,url_start,url_end",
    [(Time('1995/06/03 1:00', '1995/06/05'),
      'https://umbra.nascom.nasa.gov/goes/fits/1995/go07950603.fits',
      'https://umbra.nascom.nasa.gov/goes/fits/1995/go07950605.fits'),
     (Time('2008/06/02 12:00', '2008/06/04'),
      'https://umbra.nascom.nasa.gov/goes/fits/2008/go1020080602.fits',
      'https://umbra.nascom.nasa.gov/goes/fits/2008/go1020080604.fits'),
     (Time('2020/08/02', '2020/08/04'),
      'https://data.ngdc.noaa.gov/platforms/solar-space-observing-satellites/goes/'
      'goes16/l2/data/xrsf-l2-flx1s_science/2020/08/sci_xrsf-l2-flx1s_g16_d20200802_v2-0-1.nc',
      'https://data.ngdc.noaa.gov/platforms/solar-space-observing-satellites/goes/'
      'goes17/l2/data/xrsf-l2-flx1s_science/2020/08/sci_xrsf-l2-flx1s_g17_d20200804_v2-0-1.nc')])
def test_get_url_for_time_range(LCClient, timerange, url_start, url_end):
    qresponse = LCClient.search(timerange)
    urls = [i['url'] for i in qresponse]
    assert isinstance(urls, list)
    assert urls[0] == url_start
    assert urls[-1] == url_end


@pytest.mark.remote_data
@pytest.mark.parametrize("timerange, url_start, url_end",
                         [(a.Time('1999/01/10 00:10', '1999/01/20'),
                           'https://umbra.nascom.nasa.gov/goes/fits/1999/go10990110.fits',
                           'https://umbra.nascom.nasa.gov/goes/fits/1999/go1019990120.fits')])
def test_get_overlap_urls(LCClient, timerange, url_start, url_end):
    qresponse = LCClient.search(timerange, a.goes.SatelliteNumber.ten)
    urls = [i['url'] for i in qresponse]
    assert len(urls) == 14
    assert urls[0] == url_start
    assert urls[-1] == url_end


@pytest.mark.remote_data
@pytest.mark.parametrize("timerange, url_start, url_end",
                         [(a.Time("2009/08/30 00:10", "2009/09/02"),
                           "https://umbra.nascom.nasa.gov/goes/fits/2009/go1020090830.fits",
                           "https://satdat.ngdc.noaa.gov/sem/goes/data/science/xrs/goes14/gxrs-l2-irrad_science/"
                           "2009/09/sci_gxrs-l2-irrad_g14_d20090902_v0-0-0.nc")])
def test_get_overlap_providers(LCClient, timerange, url_start, url_end):
    qresponse = LCClient.search(timerange)
    urls = [i['url'] for i in qresponse]
    assert len(urls) == 4
    assert urls[0] == url_start
    assert urls[-1] == url_end


@pytest.mark.remote_data
@pytest.mark.parametrize("timerange, url_old, url_new",
                         [(a.Time('2013/10/28', '2013/10/29'),
                           "https://umbra.nascom.nasa.gov/goes/fits/2013/go1520131028.fits",
                           "https://satdat.ngdc.noaa.gov/sem/goes/data/science/xrs/goes13/gxrs-l2-irrad_science/"
                           "2013/10/sci_gxrs-l2-irrad_g13_d20131028_v0-0-0.nc")])
def test_old_data_access(timerange, url_old, url_new):
    # test first for old data
    qr = Fido.search(timerange, a.Instrument("XRS"), a.Provider("SDAC"))
    urls = [r['url'] for r in qr.get_response(0)]
    assert urls[0] == url_old

    # now test for new data
    qr = Fido.search(timerange, a.Instrument("XRS"))
    urls = [r['url'] for r in qr.get_response(0)]
    assert urls[0] == url_new


@given(goes_time())
def test_can_handle_query(time):
    ans1 = goes.XRSClient._can_handle_query(time, Instrument('XRS'))
    assert ans1 is True
    ans2 = goes.XRSClient._can_handle_query(time)
    assert ans2 is False
    ans3 = goes.XRSClient._can_handle_query(time, Instrument('eve'))
    assert ans3 is False


@pytest.mark.filterwarnings('ignore:ERFA function.*dubious year')
@pytest.mark.remote_data
def test_fixed_satellite(LCClient):
    ans1 = LCClient.search(a.Time("2017/01/01 2:00", "2017/01/02 2:10"),
                           a.Instrument.xrs,
                           a.goes.SatelliteNumber.fifteen)

    for resp in ans1:

        assert "g15" in resp['url']

    ans1 = LCClient.search(a.Time("2017/01/01", "2017/01/02 23:00"),
                           a.Instrument.xrs,
                           a.goes.SatelliteNumber(13))

    for resp in ans1:

        assert "g13" in resp['url']

    ans1 = LCClient.search(a.Time("1999/1/13", "1999/1/16"),
                           a.Instrument.xrs,
                           a.goes.SatelliteNumber(8))

    for resp in ans1:
        assert "go08" in resp['url']


@pytest.mark.parametrize("time", [
    Time('2005/4/27', '2005/4/27 12:00'),
    Time('2016/2/4', '2016/2/10')])
@pytest.mark.remote_data
def test_query(LCClient, time):
    qr1 = LCClient.search(time, Instrument('XRS'))
    assert isinstance(qr1, QueryResponse)
    # We only compare dates here as the start time of the qr will always be the
    # start of the day.
    assert qr1.time_range().start.strftime('%Y-%m-%d') == time.start.strftime('%Y-%m-%d')

    almost_day = TimeDelta(1*u.day - 1*u.millisecond)
    end = parse_time(time.end.strftime('%Y-%m-%d')) + almost_day
    assert is_time_equal(qr1.time_range().end, end)


@pytest.mark.remote_data
@pytest.mark.parametrize("time, instrument", [
    (Time('1983/06/17', '1983/06/18'), Instrument('XRS')),
    (Time('2012/10/4', '2012/10/6'), Instrument('XRS')),
])
def test_get(LCClient, time, instrument):
    qr1 = LCClient.search(time, instrument)
    download_list = LCClient.fetch(qr1)
    assert len(download_list) == len(qr1)


@pytest.mark.remote_data
def test_new_logic(LCClient):
    qr = LCClient.search(Time('2012/10/4 20:20', '2012/10/6'), Instrument('XRS'))
    download_list = LCClient.fetch(qr)
    assert len(download_list) == len(qr)


@pytest.mark.remote_data
@pytest.mark.parametrize(
    "time, instrument, expected_num_files",
    [(a.Time("2012/10/4", "2012/10/5"), a.Instrument.goes, 4),
     (a.Time('2013-10-28 01:00', '2013-10-28 03:00'), a.Instrument('XRS'), 2)])
def test_fido(time, instrument, expected_num_files):
    qr = Fido.search(time, instrument)
    assert isinstance(qr, UnifiedResponse)
    response = Fido.fetch(qr)
    assert len(response) == qr._numfile
    assert len(response) == expected_num_files


def test_attr_reg():
    assert a.Instrument.goes == a.Instrument("GOES")
    assert a.Instrument.xrs == a.Instrument("XRS")
    assert a.goes.SatelliteNumber.two == a.goes.SatelliteNumber("2")


def test_client_repr(LCClient):
    """
    Repr check
    """
    output = str(LCClient)
    assert output[:50] == 'sunpy.net.dataretriever.sources.goes.XRSClient\n\nPr'


def mock_query_object(LCClient):
    """
    Creating a Query Response object and prefilling it with some information
    """
    # Creating a Query Response Object
    start = '2016/1/1'
    end = '2016/1/1 23:59:59'
    obj = {
        'Time': TimeRange(parse_time(start), parse_time(end)),
        'Start Time': parse_time(start),
        'End Time': parse_time(end),
        'Instrument': 'GOES',
        'Physobs': 'irradiance',
        'Source': 'GOES',
        'Provider': 'NOAA',
        'SatelliteNumber': '15',
        'url': 'https://umbra.nascom.nasa.gov/goes/fits/2016/go1520160101.fits'
    }
    results = QueryResponse([obj], client=LCClient)
    return results


def test_show(LCClient):
    mock_qr = mock_query_object(LCClient)
    qrshow0 = mock_qr.show()
    qrshow1 = mock_qr.show('Start Time', 'Instrument')
    allcols = ['Start Time', 'End Time', 'Instrument', 'Physobs', 'Source',
               'Provider', 'SatelliteNumber']
    assert qrshow0.colnames == allcols
    assert qrshow1.colnames == ['Start Time', 'Instrument']
    assert qrshow0['Instrument'][0] == 'GOES'
