# -*- coding: utf-8 -*-
"""
Created on Wed Mar 26 20:17:06 2014

@author: stuart
"""
import datetime
import astropy.time
import pytest

from sunpy.time import parse_time
from sunpy.net.jsoc import JSOCClient
from sunpy.net.vso.vso import Results

client = JSOCClient()

def test_payload():
    start = parse_time('2012/1/1T00:00:00')
    end = parse_time('2012/1/1T00:00:45')

    payload = client._make_query_payload(start, end, 'hmi.M_42s')

    payload_expected = {
       'ds':'{0}[{1}-{2}]'.format('hmi.M_42s', start.strftime("%Y.%m.%d_%H:%M:%S_TAI"),
                                       end.strftime("%Y.%m.%d_%H:%M:%S_TAI")),
       'format':'json',
       'method':'url',
       'notify':'',
       'op':'exp_request',
       'process':'n=0|no_op',
       'protocol':'FITS,compress Rice',
       'requestor':'none',
       'filenamefmt':'{0}.{{T_REC:A}}.{{CAMERA}}.{{segment}}'.format('hmi.M_42s')
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

def test_status_request():
    r = client._request_status('none')
    assert r.json() == {u'status': 4, u'error': u"Bad RequestID 'none' provided."}

@pytest.mark.online
def test_post():
    responses = client.jsoc_query('2012/1/1T00:00:00', '2012/1/1T00:00:45', 'hmi.M_45s')
    assert isinstance(responses, list)
    assert responses[0][:4] == 'JSOC'

@pytest.mark.online
def test_post_pass():
    responses = client.jsoc_query('2012/1/1T00:00:00', '2012/1/1T00:00:45', 'hmi.M_45s', return_resp=True)
    responses[0].json()['status'] == 2
    responses[0].json()['protocol'] == 'FITS,compress Rice'
    responses[0].json()['method'] == 'url'

@pytest.mark.online
def test_post_fail(recwarn):
    client.jsoc_query('2012/1/1T00:00:00', '2012/1/1T00:00:45', 'none', return_resp=True)
    w = recwarn.pop(Warning)
    assert issubclass(w.category, Warning)
    assert "Query 0 retuned status 4 with error Cannot export series 'none' - it does not exist." in str(w.message)
    assert w.filename
    assert w.lineno

@pytest.mark.online
def test_request_status_fail():
    resp = client._request_status('none')
    assert resp.json() == {u'status': 4, u'error': u"Bad RequestID 'none' provided."}

@pytest.mark.online
def test_wait_get():
    responses = client.jsoc_query('2012/1/3T00:00:00', '2012/1/3T00:00:45', 'hmi.M_45s')
    res = client.wait_get(responses[0])
    assert isinstance(res, Results)
    assert res.total == 2