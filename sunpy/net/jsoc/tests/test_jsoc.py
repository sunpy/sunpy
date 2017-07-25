# -*- coding: utf-8 -*-
"""
Created on Wed Mar 26 20:17:06 2014

@author: stuart
"""
import os
import tempfile
import datetime
import astropy.table
import astropy.time
import astropy.units as u
import pytest

from sunpy.time import parse_time
from sunpy.net.jsoc import JSOCClient, JSOCResponse
from sunpy.net.download import Results
import sunpy.net.jsoc.attrs as attrs
import sunpy.net.vso.attrs as vso_attrs

client = JSOCClient()


def test_jsocresponse_double():
    j1 = JSOCResponse(table=astropy.table.Table(data=[[1, 2, 3, 4]]))
    j1.append(astropy.table.Table(data=[[1, 2, 3, 4]]))
    assert isinstance(j1, JSOCResponse)
    assert all(j1.table == astropy.table.vstack([astropy.table.Table(
        data=[[1, 2, 3, 4]]), astropy.table.Table(data=[[1, 2, 3, 4]])]))


def test_jsocresponse_single():
    j1 = JSOCResponse(table=None)
    assert len(j1) == 0
    j1.append(astropy.table.Table(data=[[1, 2, 3, 4]]))
    assert all(j1.table == astropy.table.Table(data=[[1, 2, 3, 4]]))
    assert len(j1) == 4


def test_empty_jsoc_response():
    Jresp = JSOCResponse()
    assert Jresp.table is None
    assert Jresp.query_args is None
    assert Jresp.requestIDs is None
    assert str(Jresp) == 'None'
    assert repr(Jresp) == 'None'
    assert len(Jresp) == 0


@pytest.mark.online
def test_query():
    Jresp = client.search(
        attrs.Time('2012/1/1T00:00:00', '2012/1/1T00:01:30'),
        attrs.Series('hmi.M_45s'), vso_attrs.Sample(90 * u.second))
    assert isinstance(Jresp, JSOCResponse)
    assert len(Jresp) == 2


@pytest.mark.flaky(reruns=5)
@pytest.mark.online
def test_post_pass():
    responses = client.search(
        attrs.Time('2012/1/1T00:00:00', '2012/1/1T00:00:45'),
        attrs.Series('hmi.M_45s'), attrs.Notify('jsoc@cadair.com'))
    aa = client.request_data(responses, return_resp=True)
    tmpresp = aa._d
    assert tmpresp['protocol'] == 'fits'
    assert tmpresp['method'] == 'url'


@pytest.mark.online
def test_post_wavelength():
    responses = client.search(
        attrs.Time('2010/07/30T13:30:00', '2010/07/30T14:00:00'),
        attrs.Series('aia.lev1_euv_12s'), attrs.Wavelength(193 * u.AA) |
        attrs.Wavelength(335 * u.AA), attrs.Notify('jsoc@cadair.com'))
    aa = client.request_data(responses, return_resp=True)
    tmpresp = aa[0]._d
    print(tmpresp)
    assert tmpresp['protocol'] == 'fits'
    assert tmpresp['method'] == 'url'
    assert tmpresp['count'] == '302'
    tmpresp = aa[1]._d
    assert tmpresp['protocol'] == 'fits'
    assert tmpresp['method'] == 'url'
    assert tmpresp['count'] == '302'


@pytest.mark.online
def test_post_notify_fail():
    responses = client.search(
        attrs.Time('2012/1/1T00:00:00', '2012/1/1T00:00:45'),
        attrs.Series('hmi.M_45s'))
    with pytest.raises(ValueError):
        client.request_data(responses)


@pytest.mark.online
def test_post_wave_series():
    with pytest.raises(TypeError):
        client.search(
            attrs.Time('2012/1/1T00:00:00', '2012/1/1T00:00:45'),
            attrs.Series('hmi.M_45s') | attrs.Series('aia.lev1_euv_12s'),
            attrs.Wavelength(193 * u.AA) | attrs.Wavelength(335 * u.AA))


@pytest.mark.online
def test_wait_get():
    responses = client.search(
        attrs.Time('2012/1/1T1:00:36', '2012/1/1T01:00:38'),
        attrs.Series('hmi.M_45s'), attrs.Notify('jsoc@cadair.com'))
    path = tempfile.mkdtemp()
    res = client.fetch(responses, path=path)
    assert isinstance(res, Results)
    assert res.total == 1


@pytest.mark.online
def test_get_request():
    responses = client.search(
        attrs.Time('2012/1/1T1:00:36', '2012/1/1T01:00:38'),
        attrs.Series('hmi.M_45s'), attrs.Notify('jsoc@cadair.com'))

    bb = client.request_data(responses)
    path = tempfile.mkdtemp()
    aa = client.get_request(bb, path=path)
    assert isinstance(aa, Results)


@pytest.mark.online
def test_results_filenames():
    responses = client.search(
        attrs.Time('2014/1/1T1:00:36', '2014/1/1T01:01:38'),
        attrs.Series('hmi.M_45s'), attrs.Notify('jsoc@cadair.com'))
    path = tempfile.mkdtemp()
    aa = client.fetch(responses, path=path)
    assert isinstance(aa, Results)
    files = aa.wait(progress=False)
    assert len(files) == len(responses)
    for hmiurl in aa.map_:
        assert os.path.isfile(hmiurl)


@pytest.mark.online
def test_invalid_query():
    with pytest.raises(ValueError):
        client.search(attrs.Time('2012/1/1T01:00:00', '2012/1/1T01:00:45'))


@pytest.mark.online
def test_make_recordset():
    d1 = {'end_time': datetime.datetime(2014, 1, 1, 1, 0, 35),
          'series': 'aia.lev1_euv_12s',
          'start_time': datetime.datetime(2014, 1, 1, 0, 0, 35)
         }
    r1 = 'aia.lev1_euv_12s[2014.01.01_00:00:35_TAI-2014.01.01_01:00:35_TAI]'
    assert client._make_recordset(**d1) == r1
    d1.update({'segment': 'image'})
    r2 = 'aia.lev1_euv_12s[2014.01.01_00:00:35_TAI-2014.01.01_01:00:35_TAI]{image}'
    assert client._make_recordset(**d1) == r2
    d1.update({'segment': ['image'], 'wavelength': 304*u.AA})
    r3 = 'aia.lev1_euv_12s[2014.01.01_00:00:35_TAI-2014.01.01_01:00:35_TAI][304]{image}'
    assert client._make_recordset(**d1) == r3

    d2 = {'end_time': datetime.datetime(2014, 1, 1, 1, 0, 35),
          'series': 'hmi.sharp_720s',
          'start_time': datetime.datetime(2014, 1, 1, 0, 0, 35)
         }
    d2.update({'primekey': {'HARPNUM': '4864'}})
    r4 = 'hmi.sharp_720s[4864][2014.01.01_00:00:35_TAI-2014.01.01_01:00:35_TAI]'
    assert client._make_recordset(**d2) == r4
    d2['primekey'] = {'HARPNUM': '4864','Foo': '123'}
    assert client._make_recordset(**d2) == r4


@pytest.mark.online
def test_lookup_records_errors():
    d1 = {'end_time': datetime.datetime(2014, 1, 1, 1, 0, 35),
          'start_time': datetime.datetime(2014, 1, 1, 0, 0, 35)
         }
    with pytest.raises(ValueError):
        client._lookup_records(d1)

    d1.update({'series': 'aia.lev1_euv_12s'})
    d1.update({'keys' : 123})
    with pytest.raises(ValueError):
        client._lookup_records(d1)

    d1['keys'] = 'T_OBS'
    d1.update({'segment': 123})
    with pytest.raises(ValueError):
        client._lookup_records(d1)

    d1['segment'] = 'image'
    d1.update({'primekey': {'foo': 'bar'}})
    with pytest.raises(TypeError):
        client._lookup_records(d1)

    d1.update({'segment': 'foo'})
    with pytest.raises(TypeError):
        client._lookup_records(d1)
