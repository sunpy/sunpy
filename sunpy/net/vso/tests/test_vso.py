import pytest
from parfive import Results

import astropy.units as u

from sunpy.net import _attrs as core_attrs
from sunpy.net import attr
from sunpy.net import attrs as a
from sunpy.net.vso import attrs as va
from sunpy.net.vso.legacy_response import QueryResponse
from sunpy.net.vso.table_response import VSOQueryResponseTable, iter_sort_response
from sunpy.net.vso.vso import VSOClient, build_client, get_online_vso_url
from sunpy.tests.mocks import MockObject
from sunpy.time import parse_time


class MockQRRecord:
    """
    Used to test sunpy.net.vso.QueryResponse.build_table(...)
    """
    def __new__(cls, start_time=None, end_time=None, size=0, source='SOHO', instrument='aia',
                extent=None, fileid="spam"):
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
    qr = client.search(
        core_attrs.Time('2011-06-07 06:33', '2011-06-07 06:33:08'),
        core_attrs.Instrument('aia'), core_attrs.Wavelength(171 * u.AA),
        response_format="table")
    tmp_dir = tmpdir / "{file}"
    files = client.fetch(qr, path=tmp_dir)

    assert len(files) == 1

    # The construction of a VSO filename is bonkers complex, so there is no
    # practical way to determine what it should be in this test, so we just
    # put it here.
    assert "aia_lev1_171a_2011_06_07t06_33_02_77z_image_lev1.fits" in files[0]


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
    qr = client.search(stereo, response_format="table")
    downloader = MockDownloader()
    res = client.fetch(qr, wait=False, downloader=downloader)
    assert downloader.download_called is False
    assert res == Results()


def test_non_str_instrument():
    # Sanity Check
    assert isinstance(core_attrs.Instrument("lyra"), core_attrs.Instrument)

    with pytest.raises(ValueError):
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
    When explicitly suppling an 'Extent' only the 'type' is stored
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
    with pytest.raises(ConnectionError):
        VSOClient()


def test_build_client(mocker):
    mocker.patch('sunpy.net.vso.vso.check_connection', return_value=None)
    with pytest.raises(ConnectionError):
        build_client(url="http://notathing.com/", port_name="spam")


def test_build_client_params():
    with pytest.raises(ValueError):
        build_client(url="http://notathing.com/")


@pytest.mark.remote_data
def test_incorrect_content_disposition(client):
    results = client.search(
        core_attrs.Time('2011/1/1 01:00', '2011/1/1 01:02'),
        core_attrs.Instrument('mdi'), response_format="table")
    files = client.fetch(results[0:1])

    assert len(files) == 1
    assert files[0].endswith("mdi_vw_v_9466622_9466622.tar")
    assert "Content" not in files[0]


@pytest.mark.parametrize("query, handle", [
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

    print(legacy_response)
    print(table_response)


def test_row_to_table(mocker, mock_build_client, client, mock_table_response):
    mock_table_response.client = client
    # we want to assert that as_table is being called, but if it returns an
    # empty list the rest of the fetch method (which does network stuff) is
    # skipped.
    as_table = mocker.patch("sunpy.net.base_client.QueryResponseRow.as_table", return_value=[])
    client.fetch(mock_table_response[0])
    assert as_table.called
