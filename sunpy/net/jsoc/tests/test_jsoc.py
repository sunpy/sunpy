# -*- coding: utf-8 -*-
"""
Created on Wed Mar 26 20:17:06 2014

@author: stuart
"""
import time
import datetime
import astropy.table
import astropy.time
import astropy.units as u
import pytest

from sunpy.time import parse_time
from sunpy.net.jsoc import JSOCClient, JSOCResponse
from sunpy.net.download import Results
import sunpy.net.jsoc.attrs as attrs


client = JSOCClient()


def test_jsocresponse_double():
    j1 = JSOCResponse(table=astropy.table.Table(data=[[1,2,3,4]]))
    j1.append(astropy.table.Table(data=[[1,2,3,4]]))
    assert isinstance(j1, JSOCResponse)
    assert all(j1.table == astropy.table.vstack([astropy.table.Table(data=[[1,2,3,4]]),
                                                 astropy.table.Table(data=[[1,2,3,4]])]))

def test_jsocresponse_single():
    j1 = JSOCResponse(table=None)
    assert len(j1) == 0
    j1.append(astropy.table.Table(data=[[1,2,3,4]]))
    assert all(j1.table == astropy.table.Table(data=[[1,2,3,4]]))
    assert len(j1) == 4

def test_payload():
    start = parse_time('2012/1/1T00:00:00')
    end = parse_time('2012/1/1T00:00:45')

    payload = client._make_query_payload(start, end, 'hmi.M_42s')

    payload_expected = {
       'ds': '{0}[{1}-{2}]'.format('hmi.M_42s',
                                   start.strftime("%Y.%m.%d_%H:%M:%S_TAI"),
                                   end.strftime("%Y.%m.%d_%H:%M:%S_TAI")),
       'format': 'json',
       'method': 'url',
       'notify': '',
       'op': 'exp_request',
       'process': 'n=0|no_op',
       'protocol': 'FITS,compress Rice',
       'requestor': 'none',
       'filenamefmt': '{0}.{{T_REC:A}}.{{CAMERA}}.{{segment}}'.format('hmi.M_42s')
       }

    assert payload == payload_expected

def test_payload_nocompression():
    start = parse_time('2012/1/1T00:00:00')
    end = parse_time('2012/1/1T00:00:45')

    payload = client._make_query_payload(start, end, 'hmi.M_42s', compression=None)

    payload_expected = {
       'ds':'{0}[{1}-{2}]'.format('hmi.M_42s', start.strftime("%Y.%m.%d_%H:%M:%S_TAI"),
                                       end.strftime("%Y.%m.%d_%H:%M:%S_TAI")),
       'format':'json',
       'method':'url',
       'notify':'',
       'op':'exp_request',
       'process':'n=0|no_op',
       'protocol':'FITS, **NONE**',
       'requestor':'none',
       'filenamefmt':'{0}.{{T_REC:A}}.{{CAMERA}}.{{segment}}'.format('hmi.M_42s')
       }

    assert payload == payload_expected

def test_payload_protocol():
    start = parse_time('2012/1/1T00:00:00')
    end = parse_time('2012/1/1T00:00:45')

    payload = client._make_query_payload(start, end, 'hmi.M_42s', protocol='as-is')

    payload_expected = {
       'ds':'{0}[{1}-{2}]'.format('hmi.M_42s', start.strftime("%Y.%m.%d_%H:%M:%S_TAI"),
                                       end.strftime("%Y.%m.%d_%H:%M:%S_TAI")),
       'format':'json',
       'method':'url',
       'notify':'',
       'op':'exp_request',
       'process':'n=0|no_op',
       'protocol':'as-is',
       'requestor':'none',
       'filenamefmt':'{0}.{{T_REC:A}}.{{CAMERA}}.{{segment}}'.format('hmi.M_42s')
       }

    assert payload == payload_expected

def test_process_time_string():
    start = client._process_time('2012/1/1T00:00:00')
    assert start == datetime.datetime(year=2012, month=1, day=1, second=34)

def test_process_time_datetime():
    start = client._process_time(datetime.datetime(year=2012, month=1, day=1))
    assert start == datetime.datetime(year=2012, month=1, day=1, second=34)

def test_process_time_astropy():
    start = client._process_time(astropy.time.Time('2012-01-01T00:00:00', format='isot', scale='utc'))
    assert start == datetime.datetime(year=2012, month=1, day=1, second=34)

def test_process_time_astropy_tai():
    start = client._process_time(astropy.time.Time('2012-01-01T00:00:00', format='isot', scale='tai'))
    assert start == datetime.datetime(year=2012, month=1, day=1, second=0)

@pytest.mark.online
def test_status_request():
    r = client._request_status('none')
    assert r.json() == {u'status': 4, u'error': u"Bad RequestID 'none' provided."}

def  test_empty_jsoc_response():
    Jresp = JSOCResponse()
    assert Jresp.table is None
    assert Jresp.query_args is None
    assert Jresp.requestIDs is None
    assert str(Jresp) == 'None'
    assert repr(Jresp) == 'None'
    assert len(Jresp) == 0

@pytest.mark.online
def test_query():
    Jresp = client.query(attrs.Time('2012/1/1T00:00:00', '2012/1/1T00:00:45'),
                         attrs.Series('hmi.M_45s'))
    assert isinstance(Jresp, JSOCResponse)
    assert len(Jresp) == 2


@pytest.mark.online
def test_post_pass():
    responses = client.query(attrs.Time('2012/1/1T00:00:00', '2012/1/1T00:00:45'),
                             attrs.Series('hmi.M_45s'))
    aa = client.request_data(responses, return_resp=True)
    tmpresp = aa[0].json()
    assert tmpresp['status'] == 2
    assert tmpresp['protocol'] == 'FITS,compress Rice'
    assert tmpresp['method'] == 'url'


@pytest.mark.online
def test_post_wavelength():
    responses = client.query(attrs.Time('2010/07/30T13:30:00','2010/07/30T14:00:00'),attrs.Series('aia.lev1_euv_12s'),
                             attrs.Wavelength(193*u.AA)|attrs.Wavelength(335*u.AA))
    aa = client.request_data(responses, return_resp=True)
    tmpresp = aa[0].json()
    assert tmpresp['status'] == 2
    assert tmpresp['protocol'] == 'FITS,compress Rice'
    assert tmpresp['method'] == 'url'
    assert tmpresp['rcount'] == 302

@pytest.mark.online()
def test_post_wave_series():
    with pytest.raises(TypeError):
        client.query(attrs.Time('2012/1/1T00:00:00', '2012/1/1T00:00:45'),
                     attrs.Series('hmi.M_45s')|attrs.Series('aia.lev1_euv_12s'),
                     attrs.Wavelength(193*u.AA)|attrs.Wavelength(335*u.AA))


@pytest.mark.online
def test_post_fail(recwarn):
    res = client.query(attrs.Time('2012/1/1T00:00:00', '2012/1/1T00:00:45'), attrs.Series('none'))
    client.request_data(res, return_resp=True)
    w = recwarn.pop(Warning)
    assert issubclass(w.category, Warning)
    assert "Query 0 retuned status 4 with error Cannot export series 'none' - it does not exist." in str(w.message)
    assert w.filename
    assert w.lineno

@pytest.mark.online
def test_request_status_fail():
    resp = client._request_status('none')
    assert resp.json() == {u'status': 4, u'error': u"Bad RequestID 'none' provided."}
    resp = client._request_status(['none'])
    assert resp.json() == {u'status': 4, u'error': u"Bad RequestID 'none' provided."}


@pytest.mark.online
@pytest.mark.xfail
def test_wait_get():
    responses = client.query(attrs.Time('2012/1/3T00:00:00', '2012/1/3T00:00:45'), attrs.Series( 'hmi.M_45s'))
    res = client.get(responses)
    assert isinstance(res, Results)
    assert res.total == 2

@pytest.mark.online
@pytest.mark.xfail
def test_check_request():
    responses = client.query(attrs.Time('2012/1/1T01:00:00', '2012/1/1T01:00:45'),
                             attrs.Series('hmi.M_45s'))

    bb = client.request_data(responses)
    aa = client.check_request(bb)
    assert aa == [6] or aa == [1] #Incase JSOC is being very efficient
    time.sleep(5)
    aa = client.check_request(bb)
    assert aa == [1]

@pytest.mark.online
def test_get_request():
    responses = client.query(attrs.Time('2012/1/1T01:00:00', '2012/1/1T01:00:45'),
                             attrs.Series('hmi.M_45s'))

    bb = client.request_data(responses)
    aa = client.get_request(bb)
    assert isinstance(aa, Results)

@pytest.mark.online
def test_invalid_query():
    with pytest.raises(ValueError):
        resp = client.query(attrs.Time('2012/1/1T01:00:00', '2012/1/1T01:00:45'))
