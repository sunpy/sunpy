import re
import tempfile

import pytest
from parfive import Results

import astropy.table
import astropy.time
import astropy.units as u
from astropy.coordinates import SkyCoord

import sunpy.map
import sunpy.net.attrs as a
from sunpy.net.jsoc import JSOCClient, JSOCResponse
from sunpy.util.exceptions import SunpyUserWarning

# Ensure all JSOC tests are run on the same parallel worker
pytestmark = [pytest.mark.xdist_group(name="jsoc")]


@pytest.fixture
def client():
    return JSOCClient()


@pytest.fixture
def jsoc_response_double(jsoc_test_email):
    resp = JSOCResponse([{'T_REC': '2011/01/01T00:00', 'INSTRUME': 'AIA'},
                         {'T_REC': '2011/01/02T00:00', 'INSTRUME': 'AIA'}])
    resp.query_args = [{'start_time': astropy.time.Time("2020-01-01T01:00:36.000"),
                        'end_time': astropy.time.Time("2020-01-01T01:00:38.000"),
                        'series': 'hmi.M_45s', 'notify': jsoc_test_email}]
    return resp


def test_jsocresponse_single():
    j1 = JSOCResponse(data=[[1, 2, 3, 4]])
    assert all(j1 == astropy.table.Table(data=[[1, 2, 3, 4]]))
    assert len(j1) == 4


def test_empty_jsoc_response():
    Jresp = JSOCResponse()
    assert len(Jresp) == 0
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


def test_show(client):
    jdict = {'TELESCOP': ['SDO/HMI', 'SDO/AIA'], 'CAR_ROT': [2145, 2145]}
    responses = JSOCResponse(astropy.table.Table(jdict))
    showtable = responses.show('TELESCOP')
    assert isinstance(showtable, astropy.table.Table)
    assert showtable.colnames == ['TELESCOP']
    assert showtable['TELESCOP'][0] == 'SDO/HMI'


@pytest.mark.remote_data
def test_post_notify_fail(client):
    responses = client.search(
        a.Time('2020/1/1T00:00:00', '2020/1/1T00:00:45'),
        a.jsoc.Series('hmi.M_45s'))
    with pytest.raises(ValueError, match="Email address is invalid or not registered"):
        client.request_data(responses)


@pytest.mark.remote_data()
def test_post_wave_series(client):
    with pytest.raises(TypeError, match="The series hmi.M_45s does not support wavelength attribute."):
        client.search(
            a.Time('2020/1/1T00:00:00', '2020/1/1T00:00:45'),
            a.jsoc.Series('hmi.M_45s') | a.jsoc.Series('aia.lev1_euv_12s'),
            a.Wavelength(193 * u.AA) | a.Wavelength(335 * u.AA))


@pytest.mark.remote_data
def test_wait_get(client, jsoc_test_email):
    responses = client.search(
        a.Time('2020/1/1T1:00:36', '2020/1/1T01:00:38'),
        a.jsoc.Series('hmi.M_45s'), a.jsoc.Notify(jsoc_test_email))
    path = tempfile.mkdtemp()
    res = client.fetch(responses, path=path)
    assert isinstance(res, Results)
    assert len(res) == 1


@pytest.mark.remote_data
def test_get_request_tar(client, jsoc_test_email):
    responses = client.search(
        a.Time('2012/1/1T1:00:36', '2012/1/1T01:00:38'),
        a.jsoc.Series('hmi.M_45s'), a.jsoc.Notify(jsoc_test_email))
    bb = client.request_data(responses, method='url-tar')
    bb.wait()
    path = tempfile.mkdtemp()
    aa = client.get_request(bb, path=path)
    assert isinstance(aa, Results)

    responses = client.search(
        a.Time('2012/1/1T1:00:36', '2012/1/1T01:00:38'),
        a.jsoc.Series('hmi.M_45s'), a.jsoc.Notify(jsoc_test_email),
        a.jsoc.Protocol('as-is'))
    bb = client.request_data(responses, method='url-tar')
    bb.wait()
    path = tempfile.mkdtemp()
    aa = client.get_request(bb, path=path)
    assert isinstance(aa, Results)


@pytest.mark.remote_data
def test_invalid_query(client):
    with pytest.raises(ValueError, match="Series must be specified for a JSOC Query"):
        client.search(a.Time('2020/1/1T01:00:00', '2020/1/1T01:00:45'))


@pytest.mark.remote_data
def test_lookup_records_errors(client):
    d1 = {'end_time': astropy.time.Time('2020-01-01 01:00:35'),
          'start_time': astropy.time.Time('2020-01-01 00:00:35')}
    with pytest.raises(ValueError, match="Series must be specified for a JSOC Query"):
        client._lookup_records(d1)

    d1['series'] = 'aia.lev1_euv_12s'
    d1['keys'] = 123
    with pytest.raises(TypeError, match="Keywords can only be passed as a list or comma-separated strings."):
        client._lookup_records(d1)
    d1['keys'] = 'T_OBS'
    d1['primekey'] = {'foo': 'bar'}
    with pytest.raises(ValueError, match=re.escape("Unexpected PrimeKeys were passed. The series aia.lev1_euv_12s supports the following Keywords: ['T_REC', 'WAVELNTH']")):
        client._lookup_records(d1)

    del d1['primekey']
    d1['segment'] = 123
    d1['wavelength'] = 304*u.AA
    with pytest.raises(TypeError, match="Segments can only be passed as a comma-separated string or a list of strings."):
        client._lookup_records(d1)

    d1['segment'] = 'foo'
    with pytest.raises(ValueError, match=re.escape("Unexpected Segments were passed. The series aia.lev1_euv_12s contains the following Segments ['image', 'spikes']")):
        client._lookup_records(d1)

    del d1['segment']
    d1['series'] = 'hmi.m_45s'
    # The series does not support wavelength attribute.
    with pytest.raises(TypeError, match=re.escape("The series hmi.m_45s does not support wavelength attribute. The following primekeys are supported ['T_REC', 'CAMERA']")):
        client._lookup_records(d1)


@pytest.mark.remote_data
def test_make_recordset_errors(client):
    d1 = {'series': 'aia.lev1_euv_12s'}
    with pytest.raises(ValueError, match="Time, Wavelength or an explicit PrimeKey must be specified."):
        client._make_recordset(**d1)

    d1 |= {
        'end_time': astropy.time.Time('2020-01-01 01:00:35', scale='tai'),
        'start_time': astropy.time.Time('2020-01-01 00:00:35', scale='tai'),
        'primekey': {
            'T_REC': '2020.01.01_00:00:35_TAI-2020.01.01_01:00:35_TAI'
        },
    }

    with pytest.raises(ValueError, match=re.escape("Time attribute has been passed both as a Time() and PrimeKey(). Please provide any one of them or separate them by OR operator.")):
        client._make_recordset(**d1)

    d1 |= {
        'end_time': astropy.time.Time('2020-01-01 01:00:35', scale='tai'),
        'start_time': astropy.time.Time('2020-01-01 00:00:35', scale='tai'),
        'wavelength': 604*u.AA,
        'primekey': {'WAVELNTH': '604'},
    }

    with pytest.raises(ValueError, match=re.escape("Wavelength attribute has been passed both as a Wavelength() and PrimeKey(). Please provide any one of them or separate them by OR operator.")):
        client._make_recordset(**d1)


@pytest.mark.remote_data
def test_make_recordset(client):
    d1 = {'series': 'aia.lev1_euv_12s',
          'end_time': astropy.time.Time('2020-01-01 01:00:35', scale='tai'),
          'start_time': astropy.time.Time('2020-01-01 00:00:35', scale='tai')
          }
    exp = 'aia.lev1_euv_12s[2020.01.01_00:00:35_TAI-2020.01.01_01:00:35_TAI]'
    assert client._make_recordset(**d1) == exp

    d1['wavelength'] = 604*u.AA
    exp = 'aia.lev1_euv_12s[2020.01.01_00:00:35_TAI-2020.01.01_01:00:35_TAI][604]'
    assert client._make_recordset(**d1) == exp

    del d1['wavelength']
    d1['primekey'] = ({'WAVELNTH': '604'})
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
    exp = 'hmi.sharp_720s[4864][2020.01.01_00:00:35_TAI-2020.01.01_01:00:35_TAI]{continuum,magnetogram}'
    assert client._make_recordset(**d1) == exp

    d1['sample'] = 300.0
    exp = 'hmi.sharp_720s[][2020.01.01_00:00:35_TAI-2020.01.01_01:00:35_TAI@300.0s]{continuum,magnetogram}'
    assert client._make_recordset(**d1) == exp


@pytest.mark.remote_data
def test_request_data_error(client, jsoc_test_email):
    responses = client.search(
        a.Time('2020/1/1T1:00:36', '2020/1/1T01:00:38'),
        a.jsoc.Series('hmi.M_45s'), a.jsoc.Notify(jsoc_test_email),
        a.jsoc.Protocol('foo'))
    with pytest.raises(TypeError, match="Protocols other than as-is,fits are not supported."):
        client.request_data(responses)


@pytest.mark.remote_data
def test_request_data_method(client, jsoc_test_email):
    responses = client.search(
        a.Time('2012/1/1T1:00:36', '2012/1/1T01:00:38'),
        a.jsoc.Series('hmi.M_45s'), a.jsoc.Notify(jsoc_test_email))
    req = client.request_data(responses, method='url-tar')
    req.wait()
    assert req.status == 0
    assert req._d['method'] == 'url-tar'
    assert req._d['protocol'] == 'fits'

    responses = client.search(
        a.Time('2012/1/1T1:00:36', '2012/1/1T01:00:38'),
        a.jsoc.Series('hmi.M_45s'), a.jsoc.Notify(jsoc_test_email),
        a.jsoc.Protocol('as-is'))
    req = client.request_data(responses, method='url-tar')
    req.wait()
    assert req.status == 0
    assert req._d['method'] == 'url-tar'
    assert req._d['protocol'] == 'as-is'


def test_can_handle_query_no_series(client):
    assert not client._can_handle_query(a.Time("2020/01/02", "2020/01/03"))
    assert not client._can_handle_query(a.Wavelength(17.1*u.nm))
    assert client._can_handle_query(a.jsoc.Series("hmi.M_45s"))


def test_jsoc_attrs(client):
    attrs = client.load_jsoc_values()
    assert a.jsoc.Series in attrs.keys()
    assert a.jsoc.Segment in attrs.keys()
    assert len(attrs[a.jsoc.Series]) != 0
    assert len(attrs[a.jsoc.Segment]) != 0


@pytest.mark.flaky(reruns_delay=30)
@pytest.mark.remote_data
def test_jsoc_cutout_attrs(client, jsoc_test_email, aia171_test_map):
    cutout = a.jsoc.Cutout(
        SkyCoord(-500*u.arcsec, -275*u.arcsec, frame=aia171_test_map.coordinate_frame),
        top_right=SkyCoord(150*u.arcsec, 375*u.arcsec, frame=aia171_test_map.coordinate_frame),
        tracking=True
    )
    q = client.search(
        a.Time(aia171_test_map.date, aia171_test_map.date + 1 * u.min),
        a.Wavelength(171*u.angstrom),
        a.jsoc.Series.aia_lev1_euv_12s,
        a.jsoc.Notify(jsoc_test_email),
        a.jsoc.Segment.image,
        cutout,
    )
    req = client.request_data(q, method='url', protocol='fits')
    req.wait()
    assert req.status == 0
    files = client.get_request(req)
    assert len(files) == 6
    m = sunpy.map.Map(files, sequence=True)
    assert m.data.shape == (1085, 1085, 6)


def test_row_and_warning(mocker, client, jsoc_response_double):
    mocker.patch("sunpy.net.jsoc.jsoc.JSOCClient.get_request")
    request_data = mocker.patch("sunpy.net.jsoc.jsoc.JSOCClient.request_data")
    with pytest.warns(SunpyUserWarning):
        client.fetch(jsoc_response_double[0], sleep=0)
    request_data.assert_called_once_with(jsoc_response_double[0].as_table())


@pytest.mark.remote_data
def test_check_request_keywords(client):
    responses = client.search(
        a.Time('2020/1/1T1:00:36', '2020/1/1T01:00:38'),
        a.jsoc.Series('hmi.M_45s'), a.jsoc.Keyword("QUALITY") == 1)
    assert len(responses) == 0

    responses = client.search(
        a.Time('2020/1/1T1:00:36', '2020/1/1T01:00:38'),
        a.jsoc.Series('hmi.M_45s'), a.jsoc.Keyword("QUALITY") == 0)
    assert len(responses) == 1

    responses = client.search(
        a.Time('2020/1/1T1:00:36', '2020/1/1T01:00:38'),
        a.jsoc.Series('hmi.M_45s'), a.jsoc.Keyword("QUALITY") < 2)
    assert len(responses) == 1

    with pytest.raises(ValueError, match="Keyword: 'EXPTIME' is not supported by series:"):
        client.search(
            a.Time('2020/1/1T1:00:36', '2020/1/1T01:00:38'),
            a.jsoc.Series('hmi.M_45s'), a.jsoc.Keyword("EXPTIME") < 2)

    responses = client.search(
        a.Time('2020/1/1T1:00:36', '2020/1/1T01:00:38'),
        a.jsoc.Series('aia.lev1_euv_12s'), a.jsoc.Keyword("QUALITY") < 2, a.jsoc.Keyword("EXPTIME") < 2)
    assert len(responses) == 0

    responses = client.search(
        a.Time('2020/1/1T1:00:36', '2020/1/1T01:00:38'),
        a.jsoc.Series('aia.lev1_euv_12s'), a.jsoc.Keyword("QUALITY") < 2, a.jsoc.Keyword("EXPTIME") > 2)
    assert len(responses) == 7


def test_empty_response_fetch(client):
    # Check that results from downloading an empty response are empty
    #
    # Could do a search like below, but to save a JSOC query in the tests
    # just test against an empty response.
    #
    # response = client.search(
    #     a.Time('1990-12-12T00:00:00', '1990-12-12T0:01:00'),
    #    a.jsoc.Series.hmi_b_720s
    # )
    response = JSOCResponse()
    assert len(response) == 0
    result = client.fetch(response)
    assert len(result) == 0


@pytest.mark.remote_data
def test_string_keyword_is_converted_correctly(client):
    # Bug reported here:
    # https://community.openastronomy.org/t/how-to-select-only-fuv-iris-level-1-data-using-jsoc/900
    responses = client.search(
        a.Time("2014-01-01T00:00", "2014-01-01T00:05"),
        a.jsoc.Series("iris.lev1"), a.jsoc.Keyword("INSTRUME") == "FUV"
    )
    assert len(responses) == 6

@pytest.mark.remote_data
def test_segments_query(client):
    res = client.search(a.jsoc.PrimeKey('HARPNUM', 3604),
                        a.jsoc.Series('hmi.sharp_cea_720s'))
    seg_res = client.search(a.jsoc.PrimeKey('HARPNUM', 3604),
                        a.jsoc.Series('hmi.sharp_cea_720s'),
                        a.jsoc.Segment('Bp') & a.jsoc.Segment('magnetogram'))
    # Check that the correct segments are returned and are not returned in the original query
    assert len(res) == len(seg_res) == 1223
    assert len(res.columns) != len(seg_res.columns)
    assert len(seg_res.columns) - len(res.columns) == 2
    assert "magnetogram" not in res.keys()
    assert "Bp" not in res.keys()
    assert "magnetogram" in seg_res.keys()
    assert "Bp" in seg_res.keys()
    # SUMS directories are like /SUM93/D1044682184/S00000/Bp.fits
    # So checking the start and end is the easiest way to check the segments
    assert seg_res["magnetogram"][0].startswith("/SUM")
    assert seg_res["Bp"][0].startswith("/SUM")
    assert seg_res["magnetogram"][0].endswith("magnetogram.fits")
    assert seg_res["Bp"][0].endswith("Bp.fits")
