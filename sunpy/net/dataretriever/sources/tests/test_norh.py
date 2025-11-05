from unittest import mock

import pytest
from hypothesis import given

import astropy.units as u

import sunpy.net.dataretriever.sources.norh as norh
from sunpy.net import Fido
from sunpy.net import attrs as a
from sunpy.net.fido_factory import UnifiedResponse
from sunpy.net.tests.helpers import mock_query_object
from sunpy.net.tests.strategies import time_attr


@pytest.fixture
def LCClient():
    return norh.NoRHClient()

@pytest.mark.parametrize(("timerange", "url_start", "url_end"), [
    (a.Time('2012/4/21', '2012/4/21'),
     'https://solar.nro.nao.ac.jp/norh/data/tcx/2012/04/tca120421',
     'https://solar.nro.nao.ac.jp/norh/data/tcx/2012/04/tca120421'
     ),
    (a.Time('2012/12/1', '2012/12/2'),
     'https://solar.nro.nao.ac.jp/norh/data/tcx/2012/12/tca121201',
     'https://solar.nro.nao.ac.jp/norh/data/tcx/2012/12/tca121202'
     ),
    (a.Time('2012/3/7', '2012/3/14'),
     'https://solar.nro.nao.ac.jp/norh/data/tcx/2012/03/tca120307',
     'https://solar.nro.nao.ac.jp/norh/data/tcx/2012/03/tca120314'
     )
])
def test_get_url_for_time_range(LCClient, timerange, url_start, url_end):
    with mock.patch('sunpy.net.dataretriever.sources.norh.NoRHClient.search',
                    return_value=mock_query_object(timerange, LCClient)):
        qresponse = LCClient.search(timerange, a.Wavelength(17*u.GHz))
        urls = [i['url'] for i in qresponse]
        assert urls[0] == url_start


@given(time_attr())
def test_can_handle_query(time):
    LCClient = norh.NoRHClient()
    ans1 = LCClient._can_handle_query(time, a.Instrument.norh)
    assert ans1 is True
    ans1 = LCClient._can_handle_query(time, a.Instrument.norh,
                                      a.Wavelength(10*u.GHz))
    assert ans1 is True
    ans2 = LCClient._can_handle_query(time)
    assert ans2 is False


@pytest.mark.remote_data
@pytest.mark.parametrize(("time", "wave"), [(a.Time('2007/08/13', '2007/08/14'),a.Wavelength(17*u.GHz)),(a.Time('2007/08/13', '2007/08/14'), a.Wavelength(34*u.GHz))])
def test_query(time, wave):
    LCClient = norh.NoRHClient()
    qr1 = LCClient.search(time, a.Instrument.norh, wave)
    assert qr1[0]['Start Time'].strftime('%Y-%m-%d') == time.start.strftime('%Y-%m-%d')
    assert qr1[-1]['End Time'].strftime('%Y-%m-%d') == time.end.strftime('%Y-%m-%d')


@pytest.mark.parametrize(("time", "instrument", "wave"), [
    (a.Time('2012/10/4', '2012/10/4'), a.Instrument.norh, a.Wavelength(17*u.GHz)),
    (a.Time('2012/10/4', '2012/10/4'), a.Instrument.norh, a.Wavelength(34*u.GHz))])
def test_get(LCClient, time, instrument, wave):
    with mock.patch('sunpy.net.dataretriever.sources.norh.NoRHClient.search',
                    return_value=mock_query_object(time, LCClient)):
        qr1 = LCClient.search(time, instrument, wave)
        with mock.patch('sunpy.net.dataretriever.sources.norh.NoRHClient.fetch',
                        return_value=mock_query_object(time, LCClient)):
            download_list = LCClient.fetch(qr1)
            assert len(download_list) == len(qr1)


@pytest.mark.parametrize(
    ("time", "instrument", "wave"),
    [(a.Time('2012/10/4', '2012/10/4'), a.Instrument.norh, a.Wavelength(17*u.GHz) | a.Wavelength(34*u.GHz))])
def test_fido(tmp_path, time, instrument, wave):
    with mock.patch('sunpy.net.Fido.search',
                    return_value=UnifiedResponse(mock_query_object(time, LCClient))):
        path = tmp_path / "sub"
        qr = Fido.search(time, instrument, wave)
        assert isinstance(qr, UnifiedResponse)
        with mock.patch('sunpy.net.Fido.fetch',
                        return_value=UnifiedResponse(mock_query_object(time, LCClient))):
            response = Fido.fetch(qr, path=path)
            assert len(response) == len(qr)


def test_attr_reg():
    assert a.Instrument.norh == a.Instrument('NORH')


def test_client_repr(LCClient):
    """
    Repr check
    """
    output = str(LCClient)
    assert output[:50] == 'sunpy.net.dataretriever.sources.norh.NoRHClient\n\nP'


def test_show():
    mock_qr = mock_query_object(a.Time('2016/1/1', '2016/1/1 23:59:59'), norh.NoRHClient())
    qrshow0 = mock_qr.show()
    qrshow1 = mock_qr.show('Wavelength', 'Instrument')
    allcols = {'Start Time', 'End Time', 'Instrument', 'Source',
               'Provider', 'Wavelength', 'url'}
    assert not allcols.difference(qrshow0.colnames)
    assert qrshow1.colnames == ['Wavelength', 'Instrument']
    assert qrshow0['Instrument'][0] == 'NORH'
    assert qrshow1['Wavelength'][0] == 17*u.GHz
