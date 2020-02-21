# -*- coding: utf-8 -*-
import os
import tempfile
import pandas as pd
import astropy.table
import astropy.time
import astropy.units as u
import pytest
from parfive import Results

from sunpy.net.jsoc import JSOCClient, JSOCResponse
import sunpy.net.attrs as a

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
    assert Jresp.requests is None
    assert str(Jresp) == 'None'
    assert repr(Jresp) == 'None'
    assert len(Jresp) == 0


@pytest.mark.remote_data
def test_query():
    Jresp = client.search(
        a.Time('2012/1/1T00:00:00', '2012/1/1T00:01:30'),
        a.jsoc.Series('hmi.M_45s'), a.Sample(90 * u.second))
    assert isinstance(Jresp, JSOCResponse)
    assert len(Jresp) == 2


@pytest.mark.remote_data
def test_post_pass():
    responses = client.search(
        a.Time('2012/1/1T00:00:00', '2012/1/1T00:00:45'),
        a.jsoc.Series('hmi.M_45s'), a.jsoc.Notify('jsoc@cadair.com'))
    aa = client.request_data(responses)
    tmpresp = aa._d
    assert tmpresp['protocol'] == 'fits'
    assert tmpresp['method'] == 'url'


@pytest.mark.remote_data
def test_build_table():
    responses = client.search(
        a.Time('2012/1/1T00:00:00', '2012/1/1T00:00:45'),
        a.jsoc.Series('hmi.M_45s'), a.jsoc.Notify('jsoc@cadair.com'))
    table = responses.build_table()
    assert isinstance(table, astropy.table.Table)

    columns = ['T_REC', 'TELESCOP', 'INSTRUME', 'WAVELNTH', 'CAR_ROT']
    assert columns == table.colnames


@pytest.mark.remote_data
def test_post_wavelength():
    responses = client.search(
        a.Time('2010/07/30T13:30:00', '2010/07/30T14:00:00'),
        a.jsoc.Series('aia.lev1_euv_12s'), a.jsoc.Wavelength(193 * u.AA) |
        a.jsoc.Wavelength(335 * u.AA), a.jsoc.Notify('jsoc@cadair.com'))
    aa = client.request_data(responses)
    [r.wait() for r in aa]
    tmpresp = aa[0]._d
    assert tmpresp['protocol'] == 'fits'
    assert tmpresp['method'] == 'url'
    assert tmpresp['count'] == '302'
    tmpresp = aa[1]._d
    assert tmpresp['protocol'] == 'fits'
    assert tmpresp['method'] == 'url'
    assert tmpresp['count'] == '302'


@pytest.mark.remote_data
def test_post_notify_fail():
    responses = client.search(
        a.Time('2012/1/1T00:00:00', '2012/1/1T00:00:45'),
        a.jsoc.Series('hmi.M_45s'))
    with pytest.raises(ValueError):
        client.request_data(responses)


@pytest.mark.remote_data()
def test_post_wave_series():
    with pytest.raises(TypeError):
        client.search(
            a.Time('2012/1/1T00:00:00', '2012/1/1T00:00:45'),
            a.jsoc.Series('hmi.M_45s') | a.jsoc.Series('aia.lev1_euv_12s'),
            a.jsoc.Wavelength(193 * u.AA) | a.jsoc.Wavelength(335 * u.AA))


@pytest.mark.remote_data
def test_wait_get():
    responses = client.search(
        a.Time('2012/1/1T1:00:36', '2012/1/1T01:00:38'),
        a.jsoc.Series('hmi.M_45s'), a.jsoc.Notify('jsoc@cadair.com'))
    path = tempfile.mkdtemp()
    res = client.fetch(responses, path=path)
    assert isinstance(res, Results)
    assert len(res) == 1


@pytest.mark.remote_data
def test_get_request():
    responses = client.search(
        a.Time('2012/1/1T1:00:36', '2012/1/1T01:00:38'),
        a.jsoc.Series('hmi.M_45s'), a.jsoc.Notify('jsoc@cadair.com'))

    bb = client.request_data(responses)
    path = tempfile.mkdtemp()
    aa = client.get_request(bb, path=path)
    assert isinstance(aa, Results)


@pytest.mark.remote_data
def test_invalid_query():
    with pytest.raises(ValueError):
        client.search(a.Time('2012/1/1T01:00:00', '2012/1/1T01:00:45'))


@pytest.mark.remote_data
def test_lookup_records_errors():
    d1 = {'end_time': astropy.time.Time('2014-01-01 01:00:35'),
          'start_time': astropy.time.Time('2014-01-01 00:00:35')}
    with pytest.raises(ValueError):          # Series must be specified for a JSOC Query
        client._lookup_records(d1)

    d1.update({'series': 'aia.lev1_euv_12s'})
    d1.update({'keys': 123})
    # Keywords can only be passed as a list or comma-separated strings.
    with pytest.raises(TypeError):
        client._lookup_records(d1)

    d1['keys'] = 'T_OBS'
    d1.update({'primekey': {'foo': 'bar'}})
    with pytest.raises(ValueError):          # Unexpected PrimeKeys were passed.
        client._lookup_records(d1)

    del d1['primekey']
    d1.update({'segment': 123})
    d1.update({'wavelength': 304*u.AA})
    # Segments can only be passed as a comma-separated string or a list of strings.
    with pytest.raises(TypeError):
        client._lookup_records(d1)

    d1.update({'segment': 'foo'})
    with pytest.raises(ValueError):          # Unexpected Segments were passed.
        client._lookup_records(d1)

    del d1['segment']
    d1.update({'series': 'hmi.m_45s'})
    # The series does not support wavelength attribute.
    with pytest.raises(TypeError):
        client._lookup_records(d1)


@pytest.mark.remote_data
def test_make_recordset_errors():
    d1 = {'series': 'aia.lev1_euv_12s'}
    with pytest.raises(ValueError):
        client._make_recordset(**d1)

    d1.update({
        'end_time': astropy.time.Time('2014-01-01 01:00:35', scale='tai'),
        'start_time': astropy.time.Time('2014-01-01 00:00:35', scale='tai'),
        'primekey': {'T_REC': '2014.01.01_00:00:35_TAI-2014.01.01_01:00:35_TAI'}
    })

    with pytest.raises(ValueError):
        client._make_recordset(**d1)

    d1.update({
        'end_time': astropy.time.Time('2014-01-01 01:00:35', scale='tai'),
        'start_time': astropy.time.Time('2014-01-01 00:00:35', scale='tai'),
        'wavelength': 604*u.AA,
        'primekey': {'WAVELNTH': '604'}
    })

    with pytest.raises(ValueError):
        client._make_recordset(**d1)


@pytest.mark.remote_data
def test_make_recordset():
    d1 = {'series': 'aia.lev1_euv_12s',
          'end_time': astropy.time.Time('2014-01-01 01:00:35', scale='tai'),
          'start_time': astropy.time.Time('2014-01-01 00:00:35', scale='tai')
          }
    exp = 'aia.lev1_euv_12s[2014.01.01_00:00:35_TAI-2014.01.01_01:00:35_TAI]'
    assert client._make_recordset(**d1) == exp

    d1.update({'wavelength': 604*u.AA})
    exp = 'aia.lev1_euv_12s[2014.01.01_00:00:35_TAI-2014.01.01_01:00:35_TAI][604]'
    assert client._make_recordset(**d1) == exp

    del d1['wavelength']
    d1.update({'primekey': {'WAVELNTH': '604'}})
    assert client._make_recordset(**d1) == exp

    del d1['start_time'], d1['end_time']
    d1['primekey'].update({'T_REC': '2014.01.01_00:00:35_TAI-2014.01.01_01:00:35_TAI'})
    exp = 'aia.lev1_euv_12s[2014.01.01_00:00:35_TAI-2014.01.01_01:00:35_TAI]'
    assert client._make_recordset(**d1) == exp

    d1 = {'series': 'hmi.v_45s',
          'end_time': astropy.time.Time('2014-01-01 01:00:35', scale='tai'),
          'start_time': astropy.time.Time('2014-01-01 00:00:35', scale='tai'),
          'segment': 'foo,bar'
          }
    exp = 'hmi.v_45s[2014.01.01_00:00:35_TAI-2014.01.01_01:00:35_TAI]{foo,bar}'
    assert client._make_recordset(**d1) == exp

    d1['segment'] = ['foo', 'bar']
    assert client._make_recordset(**d1) == exp

    d1 = {'series': 'hmi.sharp_720s',
          'end_time': astropy.time.Time('2014-01-01 01:00:35', scale='tai'),
          'start_time': astropy.time.Time('2014-01-01 00:00:35', scale='tai'),
          'segment': ['continuum', 'magnetogram'],
          'primekey': {'HARPNUM': '4864'}
          }
    exp = 'hmi.sharp_720s[4864][2014.01.01_00:00:35_TAI-2014.01.01_01:00:35_TAI]'\
          '{continuum,magnetogram}'
    assert client._make_recordset(**d1) == exp

    d1.update({'sample': 300.0})
    exp = 'hmi.sharp_720s[][2014.01.01_00:00:35_TAI-2014.01.01_01:00:35_TAI@300.0s]'\
          '{continuum,magnetogram}'
    assert client._make_recordset(**d1) == exp


@pytest.mark.remote_data
def test_search_metadata():
    metadata = client.search_metadata(a.Time('2014-01-01T00:00:00', '2014-01-01T00:02:00'),
                                      a.jsoc.Series('aia.lev1_euv_12s'), a.jsoc.Wavelength(304*u.AA))
    assert isinstance(metadata, pd.DataFrame)
    assert metadata.shape == (11, 176)
    for i in metadata.index.values:
        assert (i.startswith('aia.lev1_euv_12s') and i.endswith('[304]'))


@pytest.mark.remote_data
def test_request_data_error():
    responses = client.search(
        a.Time('2012/1/1T1:00:36', '2012/1/1T01:00:38'),
        a.jsoc.Series('hmi.M_45s'), a.jsoc.Notify('jsoc@cadair.com'),
        a.jsoc.Protocol('foo'))
    with pytest.raises(TypeError):
        req = client.request_data(responses)


@pytest.mark.remote_data
def test_request_data_protocol():
    responses = client.search(
        a.Time('2012/1/1T1:00:36', '2012/1/1T01:00:38'),
        a.jsoc.Series('hmi.M_45s'), a.jsoc.Notify('jsoc@cadair.com'))
    req = client.request_data(responses)
    req.wait()
    assert req._d['method'] == 'url'
    assert req._d['protocol'] == 'fits'

    responses = client.search(
        a.Time('2012/1/1T1:00:36', '2012/1/1T01:00:38'),
        a.jsoc.Series('hmi.M_45s'), a.jsoc.Notify('jsoc@cadair.com'),
        a.jsoc.Protocol('fits'))
    req = client.request_data(responses)
    req.wait()
    assert req._d['method'] == 'url'
    assert req._d['protocol'] == 'fits'

    responses = client.search(
        a.Time('2012/1/1T1:00:36', '2012/1/1T01:00:38'),
        a.jsoc.Series('hmi.M_45s'), a.jsoc.Notify('jsoc@cadair.com'),
        a.jsoc.Protocol('as-is'))
    req = client.request_data(responses)
    req.wait()
    assert req._d['method'] == 'url_quick'
    assert req._d['protocol'] == 'as-is'


@pytest.mark.remote_data
def test_check_request():
    responses = client.search(
        a.Time('2012/1/1T1:00:36', '2012/1/1T01:00:38'),
        a.jsoc.Series('hmi.M_45s'), a.jsoc.Notify('jsoc@cadair.com'))
    req = client.request_data(responses)
    req.wait()
    assert req.status == 0


@pytest.mark.remote_data
def test_results_filenames():
    responses = client.search(
        a.Time('2014/1/1T1:00:36', '2014/1/1T01:01:38'),
        a.jsoc.Series('hmi.M_45s'), a.jsoc.Notify('jsoc@cadair.com'))
    path = tempfile.mkdtemp()
    files = client.fetch(responses, path=path)
    assert isinstance(files, Results)
    assert len(files) == len(responses)
    for hmiurl in files:
        assert os.path.isfile(hmiurl)
