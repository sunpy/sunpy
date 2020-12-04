import os
import tempfile
from unittest import mock

import pandas as pd
import pytest
from parfive import Results

import astropy.table
import astropy.time
import astropy.units as u
from astropy.coordinates import SkyCoord

import sunpy.data.test
import sunpy.map
import sunpy.net.attrs as a
from sunpy.net.jsoc import JSOCClient, JSOCResponse
from sunpy.util.exceptions import SunpyDeprecationWarning, SunpyUserWarning


@pytest.fixture
def client():
    return JSOCClient()


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
    assert len(Jresp.table) == 0
    assert Jresp.query_args is None
    assert Jresp.requests is None
    assert len(Jresp) == 0


@pytest.mark.remote_data
def test_return_query_args(client):
    res = client.search(a.jsoc.PrimeKey('HARPNUM', 3604),
                        a.jsoc.Series('hmi.sharp_cea_720s'),
                        a.jsoc.Segment('Bp') & a.jsoc.Segment('magnetogram'))
    # Because res.query_args is list that contains dict
    assert 'primekey' in res.query_args[0]


@pytest.mark.remote_data
def test_query(client):
    Jresp = client.search(
        a.Time('2020/1/1T00:00:00', '2020/1/1T00:01:30'),
        a.jsoc.Series('hmi.M_45s'), a.Sample(90 * u.second))
    assert isinstance(Jresp, JSOCResponse)
    assert len(Jresp) == 2


@pytest.mark.remote_data
def test_post_pass(client):
    responses = client.search(
        a.Time('2020/1/1T00:00:00', '2020/1/1T00:00:45'),
        a.jsoc.Series('hmi.M_45s'), a.jsoc.Notify('jsoc@cadair.com'))
    aa = client.request_data(responses)
    tmpresp = aa._d
    assert tmpresp['protocol'] == 'fits'
    assert tmpresp['method'] == 'url'


@pytest.mark.remote_data
def test_build_table(client):
    responses = client.search(
        a.Time('2020/1/1T00:00:00', '2020/1/1T00:00:45'),
        a.jsoc.Series('hmi.M_45s'), a.jsoc.Notify('jsoc@cadair.com'))
    table = responses.build_table()
    assert isinstance(table, astropy.table.Table)

    columns = ['T_REC', 'TELESCOP', 'INSTRUME', 'WAVELNTH', 'CAR_ROT']
    assert columns == table.colnames


def test_show(client):
    jdict = {'TELESCOP': ['SDO/HMI', 'SDO/AIA'], 'CAR_ROT': [2145, 2145]}
    responses = JSOCResponse(table=astropy.table.Table(jdict))
    showtable = responses.show('TELESCOP')
    assert isinstance(showtable, astropy.table.Table)
    assert showtable.colnames == ['TELESCOP']
    assert showtable['TELESCOP'][0] == 'SDO/HMI'


@pytest.mark.remote_data
def test_post_wavelength(client):
    responses = client.search(
        a.Time('2020/07/30T13:30:00', '2020/07/30T14:00:00'),
        a.jsoc.Series('aia.lev1_euv_12s'), a.Wavelength(193 * u.AA) |
        a.Wavelength(335 * u.AA), a.jsoc.Notify('jsoc@cadair.com'))
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
def test_post_notify_fail(client):
    responses = client.search(
        a.Time('2020/1/1T00:00:00', '2020/1/1T00:00:45'),
        a.jsoc.Series('hmi.M_45s'))
    with pytest.raises(ValueError):
        client.request_data(responses)


@pytest.mark.remote_data()
def test_post_wave_series(client):
    with pytest.raises(TypeError):
        client.search(
            a.Time('2020/1/1T00:00:00', '2020/1/1T00:00:45'),
            a.jsoc.Series('hmi.M_45s') | a.jsoc.Series('aia.lev1_euv_12s'),
            a.Wavelength(193 * u.AA) | a.Wavelength(335 * u.AA))


@pytest.mark.remote_data
def test_wait_get(client):
    responses = client.search(
        a.Time('2020/1/1T1:00:36', '2020/1/1T01:00:38'),
        a.jsoc.Series('hmi.M_45s'), a.jsoc.Notify('jsoc@cadair.com'))
    path = tempfile.mkdtemp()
    res = client.fetch(responses, path=path)
    assert isinstance(res, Results)
    assert len(res) == 1


@pytest.mark.remote_data
def test_get_request(client):
    responses = client.search(
        a.Time('2020/1/1T1:00:36', '2020/1/1T01:00:38'),
        a.jsoc.Series('hmi.M_45s'), a.jsoc.Notify('jsoc@cadair.com'))
    bb = client.request_data(responses)
    path = tempfile.mkdtemp()
    aa = client.get_request(bb, path=path)
    assert isinstance(aa, Results)


@pytest.mark.remote_data
def test_get_request_tar(client):
    responses = client.search(
        a.Time('2012/1/1T1:00:36', '2012/1/1T01:00:38'),
        a.jsoc.Series('hmi.M_45s'), a.jsoc.Notify('jsoc@cadair.com'))
    bb = client.request_data(responses, method='url-tar')
    bb.wait()
    path = tempfile.mkdtemp()
    aa = client.get_request(bb, path=path)
    assert isinstance(aa, Results)

    responses = client.search(
        a.Time('2012/1/1T1:00:36', '2012/1/1T01:00:38'),
        a.jsoc.Series('hmi.M_45s'), a.jsoc.Notify('jsoc@cadair.com'),
        a.jsoc.Protocol('as-is'))
    bb = client.request_data(responses, method='url-tar')
    bb.wait()
    path = tempfile.mkdtemp()
    aa = client.get_request(bb, path=path)
    assert isinstance(aa, Results)


@pytest.mark.remote_data
def test_invalid_query(client):
    with pytest.raises(ValueError):
        client.search(a.Time('2020/1/1T01:00:00', '2020/1/1T01:00:45'))


@pytest.mark.remote_data
def test_lookup_records_errors(client):
    d1 = {'end_time': astropy.time.Time('2020-01-01 01:00:35'),
          'start_time': astropy.time.Time('2020-01-01 00:00:35')}
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
def test_make_recordset_errors(client):
    d1 = {'series': 'aia.lev1_euv_12s'}
    with pytest.raises(ValueError):
        client._make_recordset(**d1)

    d1.update({
        'end_time': astropy.time.Time('2020-01-01 01:00:35', scale='tai'),
        'start_time': astropy.time.Time('2020-01-01 00:00:35', scale='tai'),
        'primekey': {'T_REC': '2020.01.01_00:00:35_TAI-2020.01.01_01:00:35_TAI'}
    })

    with pytest.raises(ValueError):
        client._make_recordset(**d1)

    d1.update({
        'end_time': astropy.time.Time('2020-01-01 01:00:35', scale='tai'),
        'start_time': astropy.time.Time('2020-01-01 00:00:35', scale='tai'),
        'wavelength': 604*u.AA,
        'primekey': {'WAVELNTH': '604'}
    })

    with pytest.raises(ValueError):
        client._make_recordset(**d1)


@pytest.mark.remote_data
def test_make_recordset(client):
    d1 = {'series': 'aia.lev1_euv_12s',
          'end_time': astropy.time.Time('2020-01-01 01:00:35', scale='tai'),
          'start_time': astropy.time.Time('2020-01-01 00:00:35', scale='tai')
          }
    exp = 'aia.lev1_euv_12s[2020.01.01_00:00:35_TAI-2020.01.01_01:00:35_TAI]'
    assert client._make_recordset(**d1) == exp

    d1.update({'wavelength': 604*u.AA})
    exp = 'aia.lev1_euv_12s[2020.01.01_00:00:35_TAI-2020.01.01_01:00:35_TAI][604]'
    assert client._make_recordset(**d1) == exp

    del d1['wavelength']
    d1.update({'primekey': {'WAVELNTH': '604'}})
    assert client._make_recordset(**d1) == exp

    del d1['start_time'], d1['end_time']
    d1['primekey'].update({'T_REC': '2020.01.01_00:00:35_TAI-2020.01.01_01:00:35_TAI'})
    exp = 'aia.lev1_euv_12s[2020.01.01_00:00:35_TAI-2020.01.01_01:00:35_TAI]'
    assert client._make_recordset(**d1) == exp

    d1 = {'series': 'hmi.v_45s',
          'end_time': astropy.time.Time('2020-01-01 01:00:35', scale='tai'),
          'start_time': astropy.time.Time('2020-01-01 00:00:35', scale='tai'),
          'segment': 'foo,bar'
          }
    exp = 'hmi.v_45s[2020.01.01_00:00:35_TAI-2020.01.01_01:00:35_TAI]{foo,bar}'
    assert client._make_recordset(**d1) == exp

    d1['segment'] = ['foo', 'bar']
    assert client._make_recordset(**d1) == exp

    d1 = {'series': 'hmi.sharp_720s',
          'end_time': astropy.time.Time('2020-01-01 01:00:35', scale='tai'),
          'start_time': astropy.time.Time('2020-01-01 00:00:35', scale='tai'),
          'segment': ['continuum', 'magnetogram'],
          'primekey': {'HARPNUM': '4864'}
          }
    exp = 'hmi.sharp_720s[4864][2020.01.01_00:00:35_TAI-2020.01.01_01:00:35_TAI]'\
          '{continuum,magnetogram}'
    assert client._make_recordset(**d1) == exp

    d1.update({'sample': 300.0})
    exp = 'hmi.sharp_720s[][2020.01.01_00:00:35_TAI-2020.01.01_01:00:35_TAI@300.0s]'\
          '{continuum,magnetogram}'
    assert client._make_recordset(**d1) == exp


@pytest.mark.remote_data
def test_search_metadata(client):
    with pytest.raises(SunpyDeprecationWarning):
        metadata = client.search_metadata(a.Time('2020-01-01T00:00:00', '2020-01-01T00:02:00'),
                                          a.jsoc.Series('aia.lev1_euv_12s'), a.Wavelength(304*u.AA))
        assert isinstance(metadata, pd.DataFrame)
        assert metadata.shape == (11, 176)
        for i in metadata.index.values:
            assert (i.startswith('aia.lev1_euv_12s') and i.endswith('[304]'))


@pytest.mark.remote_data
def test_request_data_error(client):
    responses = client.search(
        a.Time('2020/1/1T1:00:36', '2020/1/1T01:00:38'),
        a.jsoc.Series('hmi.M_45s'), a.jsoc.Notify('jsoc@cadair.com'),
        a.jsoc.Protocol('foo'))
    with pytest.raises(TypeError):
        client.request_data(responses)


@pytest.mark.remote_data
def test_request_data_protocol(client):
    responses = client.search(
        a.Time('2020/1/1T1:00:36', '2020/1/1T01:00:38'),
        a.jsoc.Series('hmi.M_45s'), a.jsoc.Notify('jsoc@cadair.com'))
    req = client.request_data(responses)
    req.wait()
    assert req._d['method'] == 'url'
    assert req._d['protocol'] == 'fits'

    responses = client.search(
        a.Time('2020/1/1T1:00:36', '2020/1/1T01:00:38'),
        a.jsoc.Series('hmi.M_45s'), a.jsoc.Notify('jsoc@cadair.com'),
        a.jsoc.Protocol('fits'))
    req = client.request_data(responses)
    req.wait()
    assert req._d['method'] == 'url'
    assert req._d['protocol'] == 'fits'

    responses = client.search(
        a.Time('2020/1/1T1:00:36', '2020/1/1T01:00:38'),
        a.jsoc.Series('hmi.M_45s'), a.jsoc.Notify('jsoc@cadair.com'),
        a.jsoc.Protocol('as-is'))
    req = client.request_data(responses)
    req.wait()
    assert req._d['method'] == 'url_quick'
    assert req._d['protocol'] == 'as-is'


@pytest.mark.remote_data
def test_request_data_method(client):
    responses = client.search(
        a.Time('2012/1/1T1:00:36', '2012/1/1T01:00:38'),
        a.jsoc.Series('hmi.M_45s'), a.jsoc.Notify('jsoc@cadair.com'))
    req = client.request_data(responses, method='url-tar')
    req.wait()
    assert req._d['method'] == 'url-tar'
    assert req._d['protocol'] == 'fits'

    responses = client.search(
        a.Time('2012/1/1T1:00:36', '2012/1/1T01:00:38'),
        a.jsoc.Series('hmi.M_45s'), a.jsoc.Notify('jsoc@cadair.com'),
        a.jsoc.Protocol('as-is'))
    req = client.request_data(responses, method='url-tar')
    req.wait()
    assert req._d['method'] == 'url-tar'
    assert req._d['protocol'] == 'as-is'


@pytest.mark.remote_data
def test_check_request(client):
    responses = client.search(
        a.Time('2020/1/1T1:00:36', '2020/1/1T01:00:38'),
        a.jsoc.Series('hmi.M_45s'), a.jsoc.Notify('jsoc@cadair.com'))
    req = client.request_data(responses)
    req.wait()
    assert req.status == 0


@pytest.mark.flaky(reruns_delay=30)
@pytest.mark.remote_data
def test_results_filenames(client):
    responses = client.search(
        a.Time('2020/1/1T1:00:36', '2020/1/1T01:01:38'),
        a.jsoc.Series('hmi.M_45s'), a.jsoc.Notify('jsoc@cadair.com'))
    path = tempfile.mkdtemp()
    files = client.fetch(responses, path=path)
    assert isinstance(files, Results)
    assert len(files) == len(responses)
    for hmiurl in files:
        assert os.path.isfile(hmiurl)


@pytest.mark.flaky(reruns_delay=30)
@pytest.mark.remote_data
def test_results_filenames_as_is(tmp_path, client):
    responses = client.search(
        a.Time('2020/1/1T1:00:36', '2020/1/1T01:01:38'),
        a.jsoc.Series('hmi.M_45s'), a.jsoc.Notify('jsoc@cadair.com'),
        a.jsoc.Protocol('as-is'))
    assert len(responses) == 2
    files = client.fetch(responses, path=tmp_path)
    assert isinstance(files, Results)
    assert len(files) == len(responses)
    for hmiurl in files:
        assert os.path.isfile(hmiurl)


def test_can_handle_query_no_series(client):
    assert not client._can_handle_query(a.Time("2020/01/02", "2020/01/03"))
    assert not client._can_handle_query(a.Wavelength(17.1*u.nm))
    assert client._can_handle_query(a.jsoc.Series("hmi.M_45s"))


@pytest.mark.remote_data
def test_max_parallel_connections(client):
    responses = client.search(
        a.Time('2020/1/1T1:00:36', '2020/1/1T01:01:38'),
        a.jsoc.Series('hmi.M_45s'), a.jsoc.Notify('jsoc@cadair.com'),
        a.jsoc.Protocol("as-is"))

    path = tempfile.mkdtemp()

    with mock.patch(
        "parfive.Downloader.download",
        new_callable=mock.MagicMock
    ) as download:

        download.side_effect = ["Mocked Downloader"]

        with pytest.warns(SunpyUserWarning):
            client.fetch(responses, path=path, max_conn=5, max_splits=5)

    assert download.called


def test_jsoc_attrs(client):
    attrs = client.load_jsoc_values()
    assert a.jsoc.Series in attrs.keys()
    assert a.jsoc.Segment in attrs.keys()
    assert len(attrs[a.jsoc.Series]) != 0
    assert len(attrs[a.jsoc.Segment]) != 0


@pytest.mark.flaky(reruns_delay=30)
@pytest.mark.remote_data
def test_jsoc_cutout_attrs(client):
    m_ref = sunpy.map.Map(sunpy.data.test.get_test_filepath('aia_171_level1.fits'))
    cutout = a.jsoc.Cutout(
        SkyCoord(-500*u.arcsec, -275*u.arcsec, frame=m_ref.coordinate_frame),
        top_right=SkyCoord(150*u.arcsec, 375*u.arcsec, frame=m_ref.coordinate_frame),
        tracking=True
    )
    q = client.search(
        a.Time(m_ref.date, m_ref.date + 1 * u.min),
        a.Wavelength(171*u.angstrom),
        a.jsoc.Series.aia_lev1_euv_12s,
        a.jsoc.Notify('jsoc@cadair.com'),  # Put your email here
        a.jsoc.Segment.image,
        cutout,
    )
    req = client.request_data(q, method='url', protocol='fits')
    req.wait()
    assert req.status == 0
    files = client.get_request(req, max_conn=2)
    assert len(files) == 6
    m = sunpy.map.Map(files, sequence=True)
    assert m.all_maps_same_shape()
    assert m.as_array().shape == (1085, 1085, 6)
