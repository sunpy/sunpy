from functools import partial
from xml.etree import ElementTree
from urllib.error import URLError, HTTPError
from urllib.request import urlopen

import numpy as np
import pytest
from parfive import Results

import astropy.units as u

from sunpy.net import _attrs as core_attrs
from sunpy.net import attr
from sunpy.net import attrs as a
from sunpy.net.vso import attrs as va
from sunpy.net.vso.legacy_response import QueryResponse
from sunpy.net.vso.table_response import VSOQueryResponseTable, iter_sort_response
from sunpy.net.vso.vso import (
    DEFAULT_URL_PORT,
    VSOClient,
    build_client,
    check_cgi_connection,
    check_connection,
    get_online_vso_url,
)
from sunpy.tests.mocks import MockObject
from sunpy.time import parse_time
from sunpy.util.exceptions import SunpyConnectionWarning, SunpyDeprecationWarning, SunpyUserWarning


@pytest.fixture(scope="session")
def check_vso_alive():
    status = [check_connection(url["url"]) for url in DEFAULT_URL_PORT]
    if not np.all(status):
        pytest.xfail("No working VSO links found.")


class MockQRRecord:
    """
    Used to test sunpy.net.vso.QueryResponse.build_table(...)
    """
    def __new__(
        cls,
        start_time=None,
        end_time=None,
        size=0,
        source='SOHO',
        instrument='aia',
        extent=None,
        fileid="spam"
    ):
        return MockObject(size=size,
                          time=MockObject(start=start_time, end=end_time),
                          source=source,
                          instrument=instrument,
                          provider="SunPy",
                          extent=extent,
                          fileid=fileid)


class MockQRResponse:
    """
    Used to test `sunpy.net.vso.vso.iter_records` and `sunpy.net.vso.vso.iter_errors`

    >>> res = MockQRResponse(items=[1, 2, 3, [4, 5]], errors=['no-connection'])  # doctest: +SKIP
    >>> res.provideritem[1].record.recorditem  # doctest: +SKIP
    [2]
    """

    def __init__(self, records=None, errors=None):

        self.provideritem = list()

        if records is not None:
            self.provideritem = [MockObject(record=MockObject(recorditem=list(records)))]

        if errors is not None:
            self.provideritem.extend([MockObject(error=err) for err in errors])


@pytest.fixture
def mock_response():
    # defining unsorted queryresult to mock test `iter_sort_response()`.
    # Incorporated cases with no None start time and without time attribute too.
    recs = [
        MockQRRecord(start_time="2021/01/01T00:00:04", fileid='t4'),
        MockQRRecord(start_time="2021/01/01T00:00:01", fileid='t1'),
        MockQRRecord(start_time="2021/01/01T00:00:02", fileid='t2'),
        MockQRRecord(start_time=None, fileid='f1'),
        MockQRRecord(start_time=None, end_time=None, fileid='f2'),
        MockQRRecord(start_time="2021/01/01T00:00:03", fileid='t3'),
    ]
    return MockQRResponse(records=recs, errors=['FAILED'])


@pytest.fixture
def mock_table_response(mock_response):
    return VSOQueryResponseTable.from_zeep_response(mock_response, client=False)


@pytest.fixture
def eit():
    return core_attrs.Instrument('eit')


@pytest.fixture
def mock_build_client(mocker):
    return mocker.patch("sunpy.net.vso.vso.build_client", return_value=True)


def test_str(mock_build_client):
    qr = VSOQueryResponseTable()
    assert str(qr) == '<No columns>'


def test_repr(mock_build_client):
    qr = VSOQueryResponseTable()
    assert '<No columns>' in repr(qr)


def test_show(mock_build_client):
    qr = VSOQueryResponseTable()
    qrshow = qr.show('Start Time', 'Source', 'Type')
    assert str(qrshow) == '<No columns>'


@pytest.mark.remote_data
def test_path(client, tmpdir):
    """
    Test that '{file}' is automatically appended to the end of a custom path if
    it is not specified.
    """
    with pytest.warns(SunpyDeprecationWarning, match="response_format"):
        qr = client.search(
            core_attrs.Time('2025-03-03 06:33', '2025-03-03 06:33:13'),
            core_attrs.Instrument('aia'), core_attrs.Wavelength(171 * u.AA),
            response_format="table")
    tmp_dir = tmpdir / "{file}"
    files = client.fetch(qr, path=tmp_dir)
    assert len(files) == 1
    # The construction of a VSO filename is BONKERS, so there is no
    # practical way to determine what it should be in this test, so we just
    # put it here.
    assert "aia.lev1.171A_2025_03_03T06_33_09.35Z.image_lev1.fits" in files[0]


@pytest.mark.filterwarnings('ignore:ERFA function.*dubious year')
@pytest.mark.remote_data
def test_no_download(client):
    """
    Test for https://github.com/sunpy/sunpy/issues/3292
    """
    class MockDownloader:
        download_called = False

        def __init__(self):
            pass

        def download(self, *args, **kwargs):
            self.download_called = True

    # this should fail
    stereo = (core_attrs.Detector('STEREO_B') &
              core_attrs.Instrument('EUVI') &
              core_attrs.Time('1900-01-01', '1900-01-01T00:10:00'))
    with pytest.warns(SunpyDeprecationWarning, match="response_format"):
        qr = client.search(stereo, response_format="table")
    downloader = MockDownloader()
    res = client.fetch(qr, wait=False, downloader=downloader)
    assert downloader.download_called is False
    assert res == Results()


def test_non_str_instrument():
    # Sanity Check
    assert isinstance(core_attrs.Instrument("lyra"), core_attrs.Instrument)
    with pytest.raises(ValueError, match="Instrument names must be strings"):
        core_attrs.Instrument(1234)


def test_iter_sort_response(mock_response):
    fileids = [i.fileid for i in iter_sort_response(mock_response)]
    # the function would have sorted records w.r.t. start time,
    # those without start time appended at last of final response.
    assert fileids == ['t1', 't2', 't3', 't4', 'f1', 'f2']


def test_from_zeep_response(mocker):
    mocker.patch("sunpy.net.vso.vso.build_client", return_value=True)
    records = (MockQRRecord(),)

    table = VSOQueryResponseTable.from_zeep_response(MockQRResponse(records), client=None)

    # These are the only None values in the table.
    source_ = table['Source']
    assert len(source_) == 1
    assert source_[0] == 'SOHO'

    instrument_ = table['Instrument']
    assert len(instrument_) == 1
    assert instrument_[0] == 'aia'

    size_ = table['Size']
    assert len(size_) == 1
    assert size_[0] == 0.0


def test_QueryResponse_build_table_with_extent_type(mocker):
    """
    When explicitly supplying an 'Extent' only the 'type' is stored
    in the built table.
    """
    mocker.patch("sunpy.net.vso.vso.build_client", return_value=True)
    e_type = MockObject(x=1.0, y=2.5, width=37, length=129.2, type='CORONA')
    table = VSOQueryResponseTable.from_zeep_response(MockQRResponse((MockQRRecord(extent=e_type),)),
                                                     client=None)
    extent = table['Extent Type'].data
    assert len(extent) == 1
    assert extent[0] == e_type.type


def test_QueryResponse_build_table_with_no_start_time(mocker):
    """
    Only the 'end' time set, no 'start' time
    """
    mocker.patch("sunpy.net.vso.vso.build_client", return_value=True)
    a_st = parse_time((2016, 2, 14, 8, 8, 12))

    records = (MockQRRecord(end_time=a_st.strftime(va._TIMEFORMAT)),)

    table = VSOQueryResponseTable.from_zeep_response(MockQRResponse(records), client=None)

    # 'End Time' is valid, there is no 'Start Time' in the table
    assert 'Start Time' not in table.columns
    end_time_ = table['End Time']
    assert len(end_time_) == 1
    assert end_time_[0].value == '2016-02-14 08:08:12.000'


def test_QueryResponse_build_table_with_no_end_time(mocker):
    """
    Only the 'start' time is set, no 'end' time
    """
    mocker.patch("sunpy.net.vso.vso.build_client", return_value=True)
    a_st = parse_time((2016, 2, 14, 8, 8, 12))

    records = (MockQRRecord(start_time=a_st.strftime(va._TIMEFORMAT)),)

    table = VSOQueryResponseTable.from_zeep_response(MockQRResponse(records), client=None)
    start_time_ = table['Start Time']
    assert len(start_time_) == 1
    assert start_time_[0].value == '2016-02-14 08:08:12.000'



@pytest.mark.remote_data
def test_vso_hmi(client, tmpdir):
    """
    This is a regression test for https://github.com/sunpy/sunpy/issues/2284
    """
    with pytest.warns(SunpyDeprecationWarning, match="response_format"):
        res = client.search(core_attrs.Time('2020-01-02 23:52:00', '2020-01-02 23:54:00'),
                            core_attrs.Instrument('HMI') | core_attrs.Instrument('AIA'), response_format="table")

    dr = client.make_getdatarequest(res)

    # Extract the DRIs from the request
    dris = dr.request.datacontainer.datarequestitem

    # 3 HMI series and one AIA
    assert len(dris) == 4

    # For each DataRequestItem assert that there is only one series in it.
    for dri in dris:
        fileids = dri.fileiditem.fileid
        series = list(map(lambda x: x.split(':')[0], fileids))
        assert all([s == series[0] for s in series])


def test_check_connection(mocker):
    mocker.patch('sunpy.net.vso.vso.urlopen',
                 side_effect=HTTPError('http://notathing.com/', 400, 'Bad Request', {}, None))
    with pytest.warns(SunpyUserWarning,
                      match='Connection to http://notathing.com/ failed with error HTTP Error 400: Bad Request.'):
        assert check_connection('http://notathing.com/') is False


def test_check_cgi_connection(mocker):
    mocker.patch('sunpy.net.vso.vso.urlopen',
                 side_effect=HTTPError('http://notathing.com/', 400, 'Bad Request', {}, None))
    with pytest.warns(SunpyUserWarning,
                      match='Connection to http://notathing.com/ failed with error HTTP Error 400: Bad Request.'):
        assert check_cgi_connection('http://notathing.com/') is False

    mocker.patch('sunpy.net.vso.vso.urlopen', side_effect=URLError('http://notathing.com/', 400))
    with pytest.warns(
            SunpyUserWarning,
            match='Connection to http://notathing.com/ failed with error <urlopen error http://notathing.com/>.'
    ):
        assert check_cgi_connection('http://notathing.com/') is False


def fail_to_open_nso_cgi(disallowed_url, url, **kwargs):
    if url == disallowed_url:
        raise URLError(disallowed_url, 404)
    return urlopen(url, **kwargs)


@pytest.mark.remote_data
def test_fallback_if_cgi_offline(check_vso_alive, mocker):  # NOQA: ARG001
    """
    This test takes the cgi endpoint URL out of the WDSL, and then disables it,
    forcing the `get_online_vso_url` function to fallback to a secondary mirror.
    """
    default_url = DEFAULT_URL_PORT[0]["url"]
    # Check that without doing anything we get the URL we expect
    mirror = get_online_vso_url()
    assert mirror["url"] == default_url

    # Open the WDSL file and extract the cgi URL
    # Doing this like this means we don't have to hard code it.
    wdsl = urlopen(default_url).read()
    t = ElementTree.fromstring(wdsl)
    ele = t.findall("{http://schemas.xmlsoap.org/wsdl/}service")[0]
    cgi_url = list(ele.iter("{http://schemas.xmlsoap.org/wsdl/soap/}address"))[0].attrib['location']

    # Now patch out that URL so we can cause it to return an error
    mocker.patch('sunpy.net.vso.vso.urlopen', side_effect=partial(fail_to_open_nso_cgi, cgi_url))

    with pytest.warns(SunpyConnectionWarning,
                      match=f"Connection to {cgi_url} failed with error .* Retrying with different url and port"):
        mirror = get_online_vso_url()
    assert mirror["url"] != "http://docs.virtualsolar.org/WSDL/VSOi_rpc_literal.wsdl"


def test_get_online_vso_url(mocker):
    """
    No wsdl links returned valid HTTP response? Return None
    """
    mocker.patch('sunpy.net.vso.vso.check_connection', return_value=None)
    assert get_online_vso_url() is None


def test_VSOClient(mocker):
    """
    Unable to find any valid VSO mirror? Raise ConnectionError
    """
    mocker.patch('sunpy.net.vso.vso.get_online_vso_url', return_value=None)
    with pytest.raises(ConnectionError, match="No online VSO mirrors could be found."):
        VSOClient()


def test_build_client(mocker):
    mocker.patch('sunpy.net.vso.vso.check_connection', return_value=None)
    with pytest.raises(ConnectionError, match="Can't connect to url http://notathing.com/"):
        build_client(url="http://notathing.com/", port_name="spam")


def test_build_client_params():
    with pytest.raises(ValueError, match="Both url and port_name must be specified if either is."):
        build_client(url="http://notathing.com/")



@pytest.mark.remote_data
@pytest.mark.filterwarnings("ignore:Can't connect to vso")
def test_incorrect_content_disposition(client):
    with pytest.warns(SunpyDeprecationWarning, match="response_format"):
        results = client.search(
            core_attrs.Time('2011/1/1 01:00', '2011/1/1 01:02'),
            core_attrs.Instrument('mdi'), response_format="table")
    files = client.fetch(results[:1])

    assert len(files) == 1
    assert files[0].endswith("mdi_vw_V_9466622_9466622.tar")
    assert "Content" not in files[0]


@pytest.mark.parametrize(('query', 'handle'), [
    ((a.Time("2011/01/01", "2011/01/02"),), True),
    ((a.Physobs.los_magnetic_field,), False),
    ((a.Time("2011/01/01", "2011/01/02"), a.Provider("SDAC"),), True),
    ((a.jsoc.Series("wibble"), a.Physobs.los_magnetic_field,), False),
])
def test_can_handle_query(query, handle):
    assert VSOClient._can_handle_query(*query) is handle


@pytest.mark.remote_data
def test_vso_attrs(client):
    """
    Check that the dict is correctly filled.
    """
    adict = client.load_vso_values()
    assert isinstance(adict, dict)
    assert len(adict.keys()) == 6
    for key, value in adict.items():
        assert isinstance(key, attr.AttrMeta)
        assert isinstance(adict[key], list)
        assert isinstance(value, list)
        for val in value:
            assert isinstance(val, list)
            assert len(val) == 2


@pytest.mark.remote_data
def test_vso_repr(client):
    """
    Repr check (it is really long)
    """
    output = str(client)
    assert output[:50] == 'sunpy.net.vso.vso.VSOClient\n\nProvides access to qu'



@pytest.mark.remote_data
def test_response_block_properties(client):
    with pytest.warns(SunpyDeprecationWarning, match="response_format"):
        res = client.search(a.Time('2020/3/4', '2020/3/6'), a.Instrument('aia'),
                            a.Wavelength(171 * u.angstrom),
                            a.Sample(10 * u.minute),
                            response_format="legacy")
    properties = res.response_block_properties()
    assert len(properties) == 0


def test_response_block_properties_table(mocker, mock_response):
    mocker.patch("sunpy.net.vso.vso.build_client", return_value=True)
    legacy_response = QueryResponse.create(mock_response)
    table_response = VSOQueryResponseTable.from_zeep_response(mock_response, client=False)

    assert str(legacy_response)
    assert str(table_response)


def test_row_to_table(mocker, mock_build_client, client, mock_table_response):
    mock_table_response.client = client
    # we want to assert that as_table is being called, but if it returns an
    # empty list the rest of the fetch method (which does network stuff) is
    # skipped.
    as_table = mocker.patch("sunpy.net.base_client.QueryResponseRow.as_table", return_value=[])
    client.fetch(mock_table_response[0])
    assert as_table.called


@pytest.mark.remote_data
def test_iris_filename(client):
    pattern = "/home/yolo/sunpy/data/{file}"
    url = "https://www.lmsal.com/solarsoft/irisa/data/level2_compressed/2018/01/02/20180102_153155_3610108077/iris_l2_20180102_153155_3610108077_SJI_1330_t000.fits.gz"
    search_results = client.search(a.Time("2018-01-02 15:31:55", "2018-01-02 15:31:55"), a.Instrument.iris)
    filename = client.mk_filename(pattern, search_results[0], None, url)
    assert filename.endswith("iris_l2_20180102_153155_3610108077_SJI_1330_t000.fits.gz")



@pytest.mark.remote_data
def test_table_noinfo_required(client):
    res = client.search(a.Time('2017/12/17 00:00:00', '2017/12/17 06:00:00'), a.Instrument('aia'), a.Wavelength(171 * u.angstrom))
    assert 'Info Required' not in res.keys()
    assert len(res) > 0


@pytest.mark.remote_data
def test_table_has_info_required_swap(client):
    res = client.search(a.Time('2020/02/15 00:00:00', '2020/02/15 20:00:00'), a.Instrument('swap'), a.Provider('ESA'), a.Source('PROBA2'))
    assert 'Info Required' in res.keys()
    assert len(res) > 0


@pytest.mark.remote_data
def test_table_has_info_required_lyra(client):
    res = client.search(a.Time('2020/02/15 00:00:00', '2020/02/17 20:00:00'), a.Instrument('lyra'), a.Provider('ESA'), a.Source('PROBA2'))
    assert 'Info Required' in res.keys()
    assert len(res) > 0


@pytest.mark.remote_data
def test_fetch_swap(client, tmp_path):
    res = client.search(a.Time('2020/02/15 00:00:00', '2020/02/15 20:00:00'), a.Instrument('swap'), a.Provider('ESA'), a.Source('PROBA2'))
    files = client.fetch(res[0:1], path=tmp_path)
    assert len(files) == 1


@pytest.mark.remote_data
def test_fetch_lyra(client, tmp_path):
    res = client.search(a.Time('2020/02/15 00:00:00', '2020/02/17 20:00:00'), a.Instrument('lyra'), a.Provider('ESA'), a.Source('PROBA2'))
    files = client.fetch(res[0:1], path=tmp_path)
    assert len(files) == 1


@pytest.mark.remote_data
def test_client_stereo_extent(client):
    res = client.search(a.Time('2008/01/14', '2008/01/14 01:00:00'), a.Instrument.secchi, a.Source('STEREO_A'), a.ExtentType('CORONA'))
    # TODO: Update when VSo Extent filtering works
    assert len(res) == 123
    assert not all(res.columns["Extent Type"] == "CORONA")


@pytest.mark.remote_data
def test_fido_stereo_extent_type(client):
    res = client.search(a.Time('2008/01/14', '2008/01/14 01:00:00'), a.Instrument.secchi, a.Source('STEREO_A'), a.ExtentType('CORONA'))
    # TODO: Update when VSo Extent filtering works
    assert len(res) == 123
    assert not all(res.columns["Extent Type"] == "CORONA")
