# -*- coding: utf-8 -*-
# Author: Florian Mayer <florian.mayer@bitsrc.org>

#pylint: disable=W0613

from __future__ import absolute_import

import tempfile
import datetime

import pytest
from six import iteritems

from astropy import units as u

from sunpy.time import TimeRange
from sunpy.net import vso
from sunpy.net.vso import attrs as va
from sunpy.net.vso import QueryResponse
from sunpy.net import attr

from sunpy.tests.mocks import MockObject


class MockQRRecord:
    """
    Used to test sunpy.net.vso.QueryResponse.build_table(...)
    """

    def __init__(self, start_time=None, end_time=None, size=0, source='SOHO', instrument='aia',
                 extent_type=None):
        self.size = size
        self.time = MockObject(start=start_time, end=end_time)
        self.source = va.Source(source)
        self.instrument = va.Instrument(instrument)
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
            self.provideritem = [MockObject(record=MockObject(recorditem=[ri])) for ri in records]

        if errors is not None:
            self.provideritem.extend([MockObject(error=err) for err in errors])


@pytest.fixture
def mock_response():
    return MockQRResponse(records=[1, 2], errors=['FAILED'])


@pytest.fixture
def eit(request):
    return va.Instrument('eit')


@pytest.fixture
def client(request):
    return vso.VSOClient()


@pytest.fixture
def iclient(request):
    return vso.InteractiveVSOClient()


def test_simpleattr_apply():
    a = attr.ValueAttr({('test', ): 1})
    dct = {}
    va.walker.apply(a, None, dct)
    assert dct['test'] == 1


def test_Time_timerange():
    t = va.Time(TimeRange('2012/1/1', '2012/1/2'))
    assert isinstance(t, va.Time)
    assert t.min == datetime.datetime(2012, 1, 1)
    assert t.max == datetime.datetime(2012, 1, 2)


def test_input_error():
    with pytest.raises(ValueError):
        va.Time('2012/1/1')


@pytest.mark.remote_data
def test_simpleattr_create(client):
    a = attr.ValueAttr({('instrument', ): 'eit'})
    assert va.walker.create(a, client.api)[0].instrument == 'eit'


def test_simpleattr_and_duplicate():
    attr = va.Instrument('foo')
    pytest.raises(TypeError, lambda: attr & va.Instrument('bar'))
    attr |= va.Source('foo')
    pytest.raises(TypeError, lambda: attr & va.Instrument('bar'))
    otherattr = va.Instrument('foo') | va.Source('foo')
    pytest.raises(TypeError, lambda: attr & otherattr)
    pytest.raises(TypeError, lambda: (attr | otherattr) & va.Instrument('bar'))
    tst = va.Instrument('foo') & va.Source('foo')
    pytest.raises(TypeError, lambda: tst & tst)


def test_simpleattr_or_eq():
    attr = va.Instrument('eit')

    assert attr | attr == attr
    assert attr | va.Instrument('eit') == attr


def test_complexattr_apply():
    tst = {('test', 'foo'): 'a', ('test', 'bar'): 'b'}
    a = attr.ValueAttr(tst)
    dct = {'test': {}}
    va.walker.apply(a, None, dct)
    assert dct['test'] == {'foo': 'a', 'bar': 'b'}


@pytest.mark.remote_data
def test_complexattr_create(client):
    a = attr.ValueAttr({('time', 'start'): 'test'})
    assert va.walker.create(a, client.api)[0].time.start == 'test'


def test_complexattr_and_duplicate():
    attr = va.Time((2011, 1, 1), (2011, 1, 1, 1))
    pytest.raises(TypeError,
                  lambda: attr & va.Time((2011, 2, 1), (2011, 2, 1, 1)))
    attr |= va.Source('foo')
    pytest.raises(TypeError,
                  lambda: attr & va.Time((2011, 2, 1), (2011, 2, 1, 1)))


def test_complexattr_or_eq():
    attr = va.Time((2011, 1, 1), (2011, 1, 1, 1))

    assert attr | attr == attr
    assert attr | va.Time((2011, 1, 1), (2011, 1, 1, 1)) == attr


def test_attror_and():
    attr = va.Instrument('foo') | va.Instrument('bar')
    one = attr & va.Source('bar')
    other = ((va.Instrument('foo') & va.Source('bar')) |
             (va.Instrument('bar') & va.Source('bar')))
    assert one == other


def test_wave_inputQuantity():
    wrong_type_mesage = "Wave inputs must be astropy Quantities"
    with pytest.raises(TypeError) as excinfo:
        va.Wavelength(10, 23)
        assert excinfo.value.message == wrong_type_mesage
    with pytest.raises(TypeError) as excinfo:
        va.Wavelength(10 * u.AA, 23)
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
        w = va.Wavelength((62 / factor) * unit, (62 / factor) * unit)
        assert int(w.min.to(u.AA, u.equivalencies.spectral()).value) == 199

    w = va.Wavelength(62 * u.eV, 62 * u.eV)
    assert int(w.min.to(u.AA, u.equivalencies.spectral()).value) == 199
    w = va.Wavelength(62e-3 * u.keV, 62e-3 * u.keV)
    assert int(w.min.to(u.AA, u.equivalencies.spectral()).value) == 199

    for factor, unit in frequency:
        w = va.Wavelength((1.506e16 / factor) * unit, (1.506e16 / factor) * unit)
        assert int(w.min.to(u.AA, u.equivalencies.spectral()).value) == 199

    w = va.Wavelength(1.506e16 * u.Hz, 1.506e16 * u.Hz)
    assert int(w.min.to(u.AA, u.equivalencies.spectral()).value) == 199
    w = va.Wavelength(1.506e7 * u.GHz, 1.506e7 * u.GHz)
    assert int(w.min.to(u.AA, u.equivalencies.spectral()).value) == 199

    with pytest.raises(u.UnitsError) as excinfo:
        va.Wavelength(10 * u.g, 23 * u.g)
    assert ('This unit is not convertable to any of [Unit("Angstrom"), Unit("kHz"), '
            'Unit("keV")]' in str(excinfo))


def test_time_xor():
    one = va.Time((2010, 1, 1), (2010, 1, 2))
    a = one ^ va.Time((2010, 1, 1, 1), (2010, 1, 1, 2))

    assert a == attr.AttrOr(
        [va.Time((2010, 1, 1), (2010, 1, 1, 1)),
         va.Time((2010, 1, 1, 2), (2010, 1, 2))])

    a ^= va.Time((2010, 1, 1, 4), (2010, 1, 1, 5))
    assert a == attr.AttrOr([
        va.Time((2010, 1, 1), (2010, 1, 1, 1)),
        va.Time((2010, 1, 1, 2), (2010, 1, 1, 4)),
        va.Time((2010, 1, 1, 5), (2010, 1, 2))
    ])


def test_wave_xor():
    one = va.Wavelength(0 * u.AA, 1000 * u.AA)
    a = one ^ va.Wavelength(200 * u.AA, 400 * u.AA)

    assert a == attr.AttrOr([va.Wavelength(0 * u.AA, 200 * u.AA),
                             va.Wavelength(400 * u.AA, 1000 * u.AA)])

    a ^= va.Wavelength(600 * u.AA, 800 * u.AA)

    assert a == attr.AttrOr(
        [va.Wavelength(0 * u.AA, 200 * u.AA), va.Wavelength(400 * u.AA, 600 * u.AA),
         va.Wavelength(800 * u.AA, 1000 * u.AA)])


def test_err_dummyattr_create():
    with pytest.raises(TypeError):
        va.walker.create(attr.DummyAttr(), None, {})


def test_err_dummyattr_apply():
    with pytest.raises(TypeError):
        va.walker.apply(attr.DummyAttr(), None, {})


def test_wave_repr():
    """Tests the __repr__ method of class vso.attrs.Wave"""
    wav = vso.attrs.Wavelength(12 * u.AA, 16 * u.AA)
    moarwav = vso.attrs.Wavelength(15 * u.AA, 12 * u.AA)
    assert repr(wav) == "<Wavelength(12.0, 16.0, 'Angstrom')>"
    assert repr(moarwav) == "<Wavelength(12.0, 15.0, 'Angstrom')>"


def test_str():
    qr = QueryResponse([])
    assert str(qr) == ('Start Time End Time Source Instrument Type\n'
                       '---------- -------- ------ ---------- ----')


def test_repr():
    qr = QueryResponse([])
    assert "Start Time End Time  Source Instrument   Type" in repr(qr)


@pytest.mark.remote_data
def test_path(client):
    """
    Test that '{file}' is automatically appended to the end of a custom path if
    it is not specified.
    """
    qr = client.search(
        va.Time('2011-06-07 06:33', '2011-06-07 06:33:08'),
        va.Instrument('aia'), va.Wavelength(171 * u.AA))
    tmp_dir = tempfile.mkdtemp()
    files = client.fetch(qr, path=tmp_dir).wait(progress=False)

    assert len(files) == 1

    # The construction of a VSO filename is bonkers complex, so there is no
    # practical way to determine what it should be in this test, so we just
    # put it here.
    assert "aia_lev1_171a_2011_06_07t06_33_02_77z_image_lev1.fits" in files[0]


def test_non_str_instrument():
    # Sanity Check
    assert isinstance(va.Instrument("lyra"), va.Instrument)

    with pytest.raises(ValueError):
        va.Instrument(1234)


@pytest.mark.parametrize("waverange, as_dict", [
    ('3 - 4 ', {'wave_wavemin': '3', 'wave_wavemax': '4', 'wave_waveunit': 'Angstrom'}),
    ('27', {'wave_wavemin': '27', 'wave_wavemax': '27', 'wave_waveunit': 'Angstrom'}),
    ('34 - 64 GHz', {'wave_wavemin': '34', 'wave_wavemax': '64', 'wave_waveunit': 'GHz'}),
    ('12-13keV', {'wave_wavemin': '12', 'wave_wavemax': '13', 'wave_waveunit': 'keV'}),
])
def test__parse_waverange(waverange, as_dict):
    assert vso.vso._parse_waverange(waverange) == as_dict


@pytest.mark.parametrize("input, expected", [
    ('12/01/2017 - 02/10/2018', dict(time_start='12/01/2017', time_end='02/10/2018')),
])
def test__parse_date(input, expected):
    assert vso.vso._parse_date(input) == expected


def test_iter_records(mock_response):
    assert list(vso.vso.iter_records(mock_response)) == [1, 2]


def test_iter_errors(mock_response):
    prov_item = list(vso.vso.iter_errors(mock_response))

    assert len(prov_item) == 1
    assert prov_item[0].error == 'FAILED'


def test_QueryResponse_build_table_defaults():
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
    assert source_[0] == str(va.Source('SOHO'))

    instrument_ = table['Instrument'].data
    assert len(instrument_) == 1
    assert instrument_[0] == str(va.Instrument('aia'))


def test_QueryResponse_build_table_with_extent_type():
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


def test_QueryResponse_build_table_with_no_start_time():
    """
    Only the 'end' time set, no 'start' time
    """
    a_st = datetime.datetime(2016, 2, 14, 8, 8, 12)

    records = (MockQRRecord(end_time=a_st.strftime(va.TIMEFORMAT)),)

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


def test_QueryResponse_build_table_with_no_end_time():
    """
    Only the 'start' time is set, no 'end' time
    """
    a_st = datetime.datetime(2016, 2, 14, 8, 8, 12)

    records = (MockQRRecord(start_time=a_st.strftime(va.TIMEFORMAT)),)

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
    res = client.search(va.Time('2017-09-02 23:52:00', '2017-09-02 23:54:00'),
                        va.Instrument('HMI') | va.Instrument('AIA'))

    dr = client.make_getdatarequest(res)

    # Extract the DRIs from the request
    dris = dr.request.datacontainer.datarequestitem

    # 3 HMI series and one AIA
    assert len(dris) == 4

    # For each DataRequestItem assert that there is only one series in it.
    for dri in dris:
        fileids = dri.fileiditem.fileid[0]
        series = list(map(lambda x: x.split(':')[0], fileids))
        assert all([s == series[0] for s in series])
