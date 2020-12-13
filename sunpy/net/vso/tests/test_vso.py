
from unittest import mock

import pytest
from parfive import Results

import astropy.units as u

from sunpy.net import _attrs as core_attrs
from sunpy.net import attr
from sunpy.net import attrs as a
from sunpy.net import vso
from sunpy.net.vso import QueryResponse
from sunpy.net.vso import attrs as va
from sunpy.net.vso.vso import HashableResponse, VSOClient, build_client, get_online_vso_url
from sunpy.tests.mocks import MockObject
from sunpy.time import TimeRange, parse_time


class MockQRRecord:
    """
    Used to test sunpy.net.vso.QueryResponse.build_table(...)
    """

    def __init__(self, start_time=None, end_time=None, size=0, source='SOHO', instrument='aia',
                 extent_type=None):
        self.size = size
        self.time = MockObject(start=start_time, end=end_time)
        self.source = a.Source(source)
        self.instrument = core_attrs.Instrument(instrument)
        self.extent = MockObject(type=None if extent_type is None else extent_type.type)


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
            self.provideritem = records

        if errors is not None:
            self.provideritem.extend([MockObject(error=err) for err in errors])


@pytest.fixture
def mock_response():
    # defining unsorted queryresult to mock test `iter_sort_response()`.
    # Incorporated cases with no None start time and without time attribute too.
    recs = [MockObject(record=MockObject(recorditem=[
        MockObject(time=MockObject(start=4), fileid='t4'),
        MockObject(time=MockObject(start=1), fileid='t1'),
        MockObject(time=MockObject(start=2), fileid='t2')
    ]))]
    rec = MockObject(record=MockObject(recorditem=[
        MockObject(time=MockObject(start=None), fileid='f1'),
        MockObject(fileid='f2'),
        MockObject(time=MockObject(start=3), fileid='t3')
    ]))
    recs.append(rec)
    return MockQRResponse(records=recs, errors=['FAILED'])


@pytest.fixture
def eit(request):
    return core_attrs.Instrument('eit')


@pytest.fixture
def client(request):
    return vso.VSOClient()


def test_simpleattr_apply():
    a = attr.ValueAttr({('test', ): 1})
    dct = {}
    va._walker.apply(a, None, dct)
    assert dct['test'] == 1


def test_Time_timerange():
    t = core_attrs.Time(TimeRange('2012/1/1', '2012/1/2'))
    assert isinstance(t, core_attrs.Time)
    assert t.min == parse_time((2012, 1, 1))
    assert t.max == parse_time((2012, 1, 2))


def test_input_error():
    with pytest.raises(ValueError):
        core_attrs.Time('2012/1/1')


@pytest.mark.remote_data
def test_simpleattr_create(client):
    a = attr.ValueAttr({('instrument', ): 'eit'})
    assert va._walker.create(a, client.api)[0].instrument == 'eit'


def test_simpleattr_and_duplicate():
    attr = core_attrs.Instrument('foo')
    pytest.raises(TypeError, lambda: attr & core_attrs.Instrument('bar'))
    attr |= a.Source('foo')
    pytest.raises(TypeError, lambda: attr & core_attrs.Instrument('bar'))
    otherattr = core_attrs.Instrument('foo') | a.Source('foo')
    pytest.raises(TypeError, lambda: attr & otherattr)
    pytest.raises(TypeError, lambda: (attr | otherattr) & core_attrs.Instrument('bar'))
    tst = core_attrs.Instrument('foo') & a.Source('foo')
    pytest.raises(TypeError, lambda: tst & tst)


def test_simpleattr_or_eq():
    attr = core_attrs.Instrument('eit')

    assert attr | attr == attr
    assert attr | core_attrs.Instrument('eit') == attr


def test_complexattr_apply():
    tst = {('test', 'foo'): 'a', ('test', 'bar'): 'b'}
    a = attr.ValueAttr(tst)
    dct = {'test': {}}
    va._walker.apply(a, None, dct)
    assert dct['test'] == {'foo': 'a', 'bar': 'b'}


@pytest.mark.remote_data
def test_complexattr_create(client):
    a = attr.ValueAttr({('time', 'start'): 'test'})
    assert va._walker.create(a, client.api)[0].time['start'] == 'test'


def test_complexattr_and_duplicate():
    attr = core_attrs.Time((2011, 1, 1), (2011, 1, 1, 1))
    pytest.raises(TypeError,
                  lambda: attr & core_attrs.Time((2011, 2, 1), (2011, 2, 1, 1)))
    attr |= a.Source('foo')
    pytest.raises(TypeError,
                  lambda: attr & core_attrs.Time((2011, 2, 1), (2011, 2, 1, 1)))


def test_complexattr_or_eq():
    attr = core_attrs.Time((2011, 1, 1), (2011, 1, 1, 1))

    assert attr | attr == attr
    assert attr | core_attrs.Time((2011, 1, 1), (2011, 1, 1, 1)) == attr


def test_attror_and():
    attr = core_attrs.Instrument('foo') | core_attrs.Instrument('bar')
    one = attr & a.Source('bar')
    other = ((core_attrs.Instrument('foo') & a.Source('bar')) |
             (core_attrs.Instrument('bar') & a.Source('bar')))
    assert one == other


def test_wave_inputQuantity():
    wrong_type_mesage = "Wave inputs must be astropy Quantities"
    with pytest.raises(TypeError) as excinfo:
        core_attrs.Wavelength(10, 23)
        assert excinfo.value.message == wrong_type_mesage
    with pytest.raises(TypeError) as excinfo:
        core_attrs.Wavelength(10 * u.AA, 23)
        assert excinfo.value.message == wrong_type_mesage


def test_wave_toangstrom():
    # TODO: this test should test that inputs are in any of spectral units
    # more than just converted to Angstroms.
    frequency = [(1, 1 * u.Hz),
                 (1e3, 1 * u.kHz),
                 (1e6, 1 * u.MHz),
                 (1e9, 1 * u.GHz)]

    energy = [(1, 1 * u.eV),
              (1e3, 1 * u.keV),
              (1e6, 1 * u.MeV)]

    for factor, unit in energy:
        w = core_attrs.Wavelength((62 / factor) * unit, (62 / factor) * unit)
        assert int(w.min.to(u.AA, u.equivalencies.spectral()).value) == 199

    w = core_attrs.Wavelength(62 * u.eV, 62 * u.eV)
    assert int(w.min.to(u.AA, u.equivalencies.spectral()).value) == 199
    w = core_attrs.Wavelength(62e-3 * u.keV, 62e-3 * u.keV)
    assert int(w.min.to(u.AA, u.equivalencies.spectral()).value) == 199

    for factor, unit in frequency:
        w = core_attrs.Wavelength((1.506e16 / factor) * unit, (1.506e16 / factor) * unit)
        assert int(w.min.to(u.AA, u.equivalencies.spectral()).value) == 199

    w = core_attrs.Wavelength(1.506e16 * u.Hz, 1.506e16 * u.Hz)
    assert int(w.min.to(u.AA, u.equivalencies.spectral()).value) == 199
    w = core_attrs.Wavelength(1.506e7 * u.GHz, 1.506e7 * u.GHz)
    assert int(w.min.to(u.AA, u.equivalencies.spectral()).value) == 199

    with pytest.raises(u.UnitsError) as excinfo:
        core_attrs.Wavelength(10 * u.g, 23 * u.g)
    assert ('This unit is not convertable to any of [Unit("Angstrom"), Unit("kHz"), '
            'Unit("keV")]' in str(excinfo.value))


def test_time_xor():
    one = core_attrs.Time((2010, 1, 1), (2010, 1, 2))
    a = one ^ core_attrs.Time((2010, 1, 1, 1), (2010, 1, 1, 2))

    assert a == attr.AttrOr(
        [core_attrs.Time((2010, 1, 1), (2010, 1, 1, 1)),
         core_attrs.Time((2010, 1, 1, 2), (2010, 1, 2))])

    a ^= core_attrs.Time((2010, 1, 1, 4), (2010, 1, 1, 5))
    assert a == attr.AttrOr([
        core_attrs.Time((2010, 1, 1), (2010, 1, 1, 1)),
        core_attrs.Time((2010, 1, 1, 2), (2010, 1, 1, 4)),
        core_attrs.Time((2010, 1, 1, 5), (2010, 1, 2))
    ])


def test_wave_xor():
    one = core_attrs.Wavelength(0 * u.AA, 1000 * u.AA)
    a = one ^ core_attrs.Wavelength(200 * u.AA, 400 * u.AA)

    assert a == attr.AttrOr([core_attrs.Wavelength(0 * u.AA, 200 * u.AA),
                             core_attrs.Wavelength(400 * u.AA, 1000 * u.AA)])

    a ^= core_attrs.Wavelength(600 * u.AA, 800 * u.AA)

    assert a == attr.AttrOr(
        [core_attrs.Wavelength(0 * u.AA, 200 * u.AA), core_attrs.Wavelength(400 * u.AA, 600 * u.AA),
         core_attrs.Wavelength(800 * u.AA, 1000 * u.AA)])


def test_err_dummyattr_create():
    with pytest.raises(TypeError):
        va._walker.create(attr.DummyAttr(), None, {})


def test_err_dummyattr_apply():
    with pytest.raises(TypeError):
        va._walker.apply(attr.DummyAttr(), None, {})


def test_wave_repr():
    """Tests the __repr__ method of class vso.attrs.Wave"""
    wav = core_attrs.Wavelength(12 * u.AA, 16 * u.AA)
    moarwav = core_attrs.Wavelength(15 * u.AA, 12 * u.AA)
    assert repr(wav) == "<sunpy.net.attrs.Wavelength(12.0, 16.0, 'Angstrom')>"
    assert repr(moarwav) == "<sunpy.net.attrs.Wavelength(12.0, 15.0, 'Angstrom')>"


@mock.patch("sunpy.net.vso.vso.build_client", return_value=True)
def test_str(mock_build_client):
    qr = QueryResponse([])
    assert str(qr) == ('Start Time End Time Source Instrument Type\n'
                       '---------- -------- ------ ---------- ----')


@mock.patch("sunpy.net.vso.vso.build_client", return_value=True)
def test_repr(mock_build_client):
    qr = QueryResponse([])
    assert "Start Time End Time Source Instrument Type" in repr(qr)


@mock.patch("sunpy.net.vso.vso.build_client", return_value=True)
def test_show(mock_build_client):
    qr = QueryResponse([])
    qrshow = qr.show('Start Time', 'Source', 'Type')
    assert str(qrshow) == ('Start Time Source Type\n'
                           '---------- ------ ----')


@pytest.mark.remote_data
def test_path(client, tmpdir):
    """
    Test that '{file}' is automatically appended to the end of a custom path if
    it is not specified.
    """
    qr = client.search(
        core_attrs.Time('2011-06-07 06:33', '2011-06-07 06:33:08'),
        core_attrs.Instrument('aia'), core_attrs.Wavelength(171 * u.AA))
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
    qr = client.search(stereo)
    downloader = MockDownloader()
    res = client.fetch(qr, wait=False, downloader=downloader)
    assert downloader.download_called is False
    assert res == Results()


def test_non_str_instrument():
    # Sanity Check
    assert isinstance(core_attrs.Instrument("lyra"), core_attrs.Instrument)

    with pytest.raises(ValueError):
        core_attrs.Instrument(1234)


@pytest.mark.parametrize("waverange, as_dict", [
    ('3 - 4 ', {'wave_wavemin': '3', 'wave_wavemax': '4', 'wave_waveunit': 'Angstrom'}),
    ('27', {'wave_wavemin': '27', 'wave_wavemax': '27', 'wave_waveunit': 'Angstrom'}),
    ('34 - 64 GHz', {'wave_wavemin': '34', 'wave_wavemax': '64', 'wave_waveunit': 'GHz'}),
    ('12-13keV', {'wave_wavemin': '12', 'wave_wavemax': '13', 'wave_waveunit': 'keV'}),
])
def test__parse_waverange(waverange, as_dict):
    assert vso.vso._parse_waverange(waverange) == as_dict


@pytest.mark.parametrize("input, expected", [
    ('12/01/2020 - 02/10/2018', dict(time_start='12/01/2020', time_end='02/10/2018')),
])
def test__parse_date(input, expected):
    assert vso.vso._parse_date(input) == expected


def test_iter_sort_response(mock_response):
    fileids = [i.fileid for i in vso.vso.iter_sort_response(mock_response)]
    # the function would have sorted records w.r.t. start time,
    # those without start time appended at last of final response.
    assert fileids == ['t1', 't2', 't3', 't4', 'f1', 'f2']


def test_iter_errors(mock_response):
    prov_item = list(vso.vso.iter_errors(mock_response))

    assert len(prov_item) == 1
    assert prov_item[0].error == 'FAILED'


@mock.patch("sunpy.net.vso.vso.build_client", return_value=True)
def test_QueryResponse_build_table_defaults(mock_build_client):
    records = (MockQRRecord(),)

    qr = vso.QueryResponse(records)
    table = qr.build_table()

    start_time_ = table['Start Time']
    assert len(start_time_) == 1
    assert start_time_[0] == 'None'

    end_time_ = table['End Time']
    assert len(end_time_) == 1
    assert end_time_[0] == 'None'

    type_ = table['Type'].data
    assert len(type_) == 1
    assert type_[0] == 'N/A'

    # Check values we did set by default in 'MockQRRecord'
    source_ = table['Source'].data
    assert len(source_) == 1
    assert source_[0][:-16] == str(a.Source('SOHO'))[:-16]

    instrument_ = table['Instrument'].data
    assert len(instrument_) == 1
    assert instrument_[0][:-16] == str(core_attrs.Instrument('aia'))[:-16]


@mock.patch("sunpy.net.vso.vso.build_client", return_value=True)
def test_QueryResponse_build_table_with_extent_type(mock_build_client):
    """
    When explcitley suppling an 'Extent' only the 'type' is stored
    in the built table.
    """
    e_type = va.Extent(x=1.0, y=2.5, width=37, length=129.2, atype='CORONA')

    qr = vso.QueryResponse((MockQRRecord(extent_type=e_type),))
    table = qr.build_table()

    extent = table['Type'].data
    assert len(extent) == 1
    assert extent[0] == e_type.type


@mock.patch("sunpy.net.vso.vso.build_client", return_value=True)
def test_QueryResponse_build_table_with_no_start_time(mock_build_client):
    """
    Only the 'end' time set, no 'start' time
    """
    a_st = parse_time((2016, 2, 14, 8, 8, 12))

    records = (MockQRRecord(end_time=a_st.strftime(va._TIMEFORMAT)),)

    qr = vso.QueryResponse(records)
    table = qr.build_table()

    start_time_ = table['Start Time']
    assert len(start_time_) == 1
    assert start_time_[0] == 'None'

    # Even though 'End Time' is valid, there is no 'Start Time'
    # marks as 'N/A'
    end_time_ = table['End Time']
    assert len(end_time_) == 1
    assert end_time_[0] == 'N/A'


@mock.patch("sunpy.net.vso.vso.build_client", return_value=True)
def test_QueryResponse_build_table_with_no_end_time(mock_build_client):
    """
    Only the 'start' time is set, no 'end' time
    """
    a_st = parse_time((2016, 2, 14, 8, 8, 12))

    records = (MockQRRecord(start_time=a_st.strftime(va._TIMEFORMAT)),)

    qr = vso.QueryResponse(records)
    table = qr.build_table()

    start_time_ = table['Start Time']
    assert len(start_time_) == 1
    assert start_time_[0] == '2016-02-14 08:08:12'

    end_time_ = table['End Time']
    assert len(end_time_) == 1
    assert end_time_[0] == 'None'


@pytest.mark.remote_data
def test_vso_hmi(client, tmpdir):
    """
    This is a regression test for https://github.com/sunpy/sunpy/issues/2284
    """
    res = client.search(core_attrs.Time('2020-01-02 23:52:00', '2020-01-02 23:54:00'),
                        core_attrs.Instrument('HMI') | core_attrs.Instrument('AIA'))

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


@mock.patch('sunpy.net.vso.vso.check_connection', return_value=None)
def test_get_online_vso_url(mock_urlopen):
    """
    No wsdl links returned valid HTTP response? Return None
    """
    assert get_online_vso_url() is None


@mock.patch('sunpy.net.vso.vso.get_online_vso_url', return_value=None)
def test_VSOClient(mock_vso_url):
    """
    Unable to find any valid VSO mirror? Raise ConnectionError
    """
    with pytest.raises(ConnectionError):
        VSOClient()


@mock.patch('sunpy.net.vso.vso.check_connection', return_value=None)
def test_build_client(mock_vso_url):
    with pytest.raises(ConnectionError):
        build_client(url="http://notathing.com/", port_name="spam")


def test_build_client_params():
    with pytest.raises(ValueError):
        build_client(url="http://notathing.com/")


@pytest.mark.remote_data
def test_vso_post_search(client):
    timerange = a.Time(('2020-01-01 00:01:05'), ('2020-01-01 00:01:10'))
    results = client.search(timerange, a.Instrument('aia') | a.Instrument('hmi'))
    trange_filter = a.Time(('2020-01-01 00:01:05'), ('2020-01-01 00:01:07'))
    wave_filter = a.Wavelength(131 * u.Angstrom)
    res_filtered = results.search(trange_filter & a.Instrument('aia') & wave_filter)
    for rec in res_filtered:
        assert parse_time(rec.time.start) <= parse_time(trange_filter.end)
        assert rec.instrument.lower() == 'aia'
        assert (float(rec.wave.wavemax) == 131.0 or float(rec.wave.wavemin) == 131.0)


@pytest.mark.remote_data
def test_incorrect_content_disposition(client):
    results = client.search(
        core_attrs.Time('2011/1/1 01:00', '2011/1/1 01:02'),
        core_attrs.Instrument('mdi'))
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
    res = client.search(a.Time('2020/3/4', '2020/3/6'), a.Instrument('aia'), a.Wavelength(171 * u.angstrom),
                        a.Sample(10 * u.minute))
    properties = res.response_block_properties()
    assert len(properties) == 0


def test_HashableResponse():
    d0 = {'instrument': 'HMI', 'wave': {'wavemin': '6173', 'wavemax': '6174'}, 'fileid': 'fid1'}
    d1 = {'instrument': 'AIA', 'wave': {'wavemin': '304', 'wavemax': '304'}, 'fileid': 'fid1'}
    d2 = {'instrument': 'AIA', 'wave': {'wavemin': '131', 'wavemax': '131'}, 'fileid': 'fid2'}
    dicts = [d0, d1, d2]
    hashRecs = list()
    for d in dicts:
        hashRecs.append(HashableResponse(d))
    assert d0['wave']['wavemax'] == hashRecs[0].wave.wavemax == '6174'
    assert d2['fileid'] == hashRecs[2].fileid
    assert hashRecs[0].__eq__(hashRecs[2]) is False
    assert hashRecs[0] == hashRecs[1]
    assert hashRecs[1].__hash__() != hashRecs[2].__hash__()
    assert len(set(hashRecs)) == 2
