import tempfile

import pytest
from hypothesis import given

import astropy.units as u

import sunpy.net.dataretriever.sources.goes as goes
from sunpy.net import Fido
from sunpy.net import attrs as a
from sunpy.net.dataretriever.client import QueryResponse
from sunpy.net.tests.strategies import time_attr
from sunpy.time import TimeRange, parse_time


@pytest.fixture
def suvi_client():
    return goes.SUVIClient()


@given(time_attr())
def test_can_handle_query(suvi_client, time):
    ans1 = suvi_client._can_handle_query(time, a.Instrument.suvi)
    assert ans1 is True
    ans2 = suvi_client._can_handle_query(time, a.Instrument.suvi,
                                         a.Wavelength(131 * u.Angstrom))
    assert ans2 is True
    ans3 = suvi_client._can_handle_query(time, a.Instrument.suvi,
                                         a.Wavelength(131 * u.Angstrom),
                                         a.Level.two)
    assert ans3 is True
    ans4 = suvi_client._can_handle_query(time)
    assert ans4 is False
    ans5 = suvi_client._can_handle_query(time, a.Instrument.aia)
    assert ans5 is False
    ans6 = suvi_client._can_handle_query(time, a.Instrument.suvi,
                                         a.goes.SatelliteNumber(16))
    assert ans6 is True


def test_get_goes_sat_num(suvi_client):
    date = parse_time('2019/06/11 00:00')
    min_satellite_number = 16  # when SUVI was first included
    assert suvi_client._get_goes_sat_num(date) >= min_satellite_number
    assert type(suvi_client._get_goes_sat_num(date)) is int


def test_get_goes_sat_num_error(suvi_client):
    date = parse_time('1800/06/11 00:00')
    with pytest.raises(ValueError):
        suvi_client._get_goes_sat_num(date)


def test_get_url_for_timerange_errors(suvi_client):
    """Check that unsupported values raise errors."""
    tr = TimeRange('2019/06/11 00:00', '2019/06/11 00:10')
    with pytest.raises(ValueError):
        suvi_client._get_url_for_timerange(tr, level=0)
    with pytest.raises(ValueError):
        suvi_client._get_url_for_timerange(tr, wavelength=100 * u.Angstrom)
    with pytest.raises(ValueError):
        suvi_client._get_url_for_timerange(tr, satellitenumber=1)


def mock_querry_object(suvi_client, start, end, wave):
    """
    Creating a Query Response object and prefilling it with some information
    """
    # Creating a Query Response Object
    obj = {
        'TimeRange': TimeRange(parse_time(start), parse_time(end)),
        'Time_start': parse_time(start),
        'Time_end': parse_time(end),
        'source': 'GOES',
        'instrument': 'SUVI',
        'physobs': 'flux',
        'provider': 'NOAA'
    }
    results = QueryResponse.create(obj, suvi_client._get_url_for_timerange(TimeRange(start, end),
                                                                           wavelength=wave), client=suvi_client)
    return results


def test_attr_reg():
    a.Instrument.suvi = a.Instrument("SUVI")
    a.goes.SatelliteNumber.A16 = a.goes.SatelliteNumber("16")


@pytest.mark.remote_data
def test_fetch_working(suvi_client):
    """
    Tests if the online server for fermi_gbm is working.
    This also checks if the mock is working well.
    """
    start = '2019/05/25 00:50'
    end = '2019/05/25 00:52'
    wave = 94 * u.Angstrom
    qr1 = suvi_client.search(a.Time(start, end), a.Instrument.suvi, a.Wavelength(wave))

    # Mock QueryResponse object
    mock_qr = mock_querry_object(suvi_client, start, end, wave)

    # Compare if two objects have the same attribute

    mock_qr = mock_qr.blocks[0]
    qr = qr1.blocks[0]

    assert mock_qr.source == qr.source
    assert mock_qr.provider == qr.provider
    assert mock_qr.physobs == qr.physobs
    assert mock_qr.instrument == qr.instrument
    assert mock_qr.url == qr.url

    assert qr1.time_range() == TimeRange("2019-05-25T00:52:00.000",
                                         "2019-05-25T00:56:00.000")

    with tempfile.TemporaryDirectory() as tmpdirname:
        download_list = suvi_client.fetch(qr1, path=tmpdirname)
    assert len(download_list) == len(qr1)


@pytest.mark.remote_data
@pytest.mark.parametrize("start, end, wave, expected_num_files",
                         [('2019/05/25 00:50', '2019/05/25 00:52', 94, 1),
                          ('2019/05/25 00:50', '2019/05/25 00:52', 131, 1),
                          ('2019/05/25 00:50', '2019/05/25 00:52', 171, 1),
                          ('2019/05/25 00:50', '2019/05/25 00:52', 195, 1),
                          ('2019/05/25 00:50', '2019/05/25 00:52', 284, 1),
                          ('2019/05/25 00:50', '2019/05/25 00:52', 304, 1)]
                         )
def test_get_url_for_time_range_level2(suvi_client, start, end, wave, expected_num_files):
    urls = suvi_client._get_url_for_timerange(TimeRange(start, end),
                                              wavelength=wave * u.Angstrom,
                                              level=2)
    assert isinstance(urls, list)
    assert len(urls) == expected_num_files


@pytest.mark.remote_data
@pytest.mark.parametrize("start, end, expected_num_files",
                         [('2019/05/25 00:50', '2019/05/25 00:52', 6)]
                         )
def test_get_url_for_time_range_level2_allwave(suvi_client, start, end, expected_num_files):
    """check that we get all wavelengths if no wavelength is given"""
    urls = suvi_client._get_url_for_timerange(TimeRange(start, end), level=2)
    assert isinstance(urls, list)
    assert len(urls) == expected_num_files


@pytest.mark.remote_data
@pytest.mark.parametrize("start, end ,wave, expected_num_files",
                         [('2019/05/25 00:50', '2019/05/25 00:54', 94, 6),
                          ('2019/05/25 00:50', '2019/05/25 00:54', 131, 3),
                          ('2019/05/25 00:50', '2019/05/25 00:54', 171, 2),
                          ('2019/05/25 00:50', '2019/05/25 00:54', 195, 7),
                          ('2019/05/25 00:50', '2019/05/25 00:54', 284, 2),
                          ('2019/05/25 00:50', '2019/05/25 00:54', 304, 4)]
                         )
def test_get_url_for_time_range_level1b(suvi_client, start, end, wave, expected_num_files):
    """check that we get all wavelengths if no wavelength is given"""
    urls = suvi_client._get_url_for_timerange(TimeRange(start, end),
                                              wavelength=wave * u.Angstrom,
                                              level='1b')
    assert isinstance(urls, list)
    assert len(urls) == expected_num_files


@pytest.mark.remote_data
@pytest.mark.parametrize("start, end ,wave, expected_num_files",
                         [('2019/05/25 00:50', '2019/05/25 00:54', 94, 6),
                          ('2019/05/25 00:50', '2019/05/25 00:54', 131, 3),
                          ('2019/05/25 00:50', '2019/05/25 00:54', 171, 2),
                          ('2019/05/25 00:50', '2019/05/25 00:54', 195, 7),
                          ('2019/05/25 00:50', '2019/05/25 00:54', 284, 2),
                          ('2019/05/25 00:50', '2019/05/25 00:54', 304, 4)]
                         )
def test_fido_onewave_level1b(start, end, wave, expected_num_files):
    result = Fido.search(a.Time(start, end), a.Instrument.suvi,
                         a.Wavelength(wave * u.Angstrom), a.Level('1b'))
    assert result.file_num == expected_num_files


@pytest.mark.remote_data
@pytest.mark.parametrize("start, end, wave1, wave2, expected_num_files",
                         [('2019/05/25 00:50', '2019/05/25 00:54', 1, 100, 6),
                          ('2019/05/25 00:50', '2019/05/25 00:54', 1, 150, 9),
                          ('2019/05/25 00:50', '2019/05/25 00:54', 1, 180, 11),
                          ('2019/05/25 00:50', '2019/05/25 00:54', 1, 200, 18),
                          ('2019/05/25 00:50', '2019/05/25 00:54', 1, 300, 20),
                          ('2019/05/25 00:50', '2019/05/25 00:54', 1, 310, 24)]
                         )
def test_fido_waverange_level1b(start, end, wave1, wave2, expected_num_files):
    """check that we get all wavelengths if no wavelength is given"""
    result = Fido.search(a.Time(start, end), a.Instrument.suvi,
                         a.Wavelength(wave1 * u.Angstrom, wave2 * u.Angstrom),
                         a.Level('1b'))
    assert result.file_num == expected_num_files


@pytest.mark.remote_data
@pytest.mark.parametrize("start, end, expected_num_files",
                         [('2019/05/25 00:50', '2019/05/25 00:52', 6)]
                         )
def test_query(suvi_client, start, end, expected_num_files):
    qr1 = suvi_client.search(a.Time(start, end), a.Instrument.suvi)
    assert isinstance(qr1, QueryResponse)
    assert len(qr1) == expected_num_files
    assert qr1.time_range().start == parse_time('2019/05/25 00:52')
    assert qr1.time_range().end == parse_time('2019/05/25 00:56')
