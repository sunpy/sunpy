import pytest
from hypothesis import given, settings

import astropy.units as u

import sunpy.net.dataretriever.sources.goes as goes
from sunpy.net import attrs as a
from sunpy.net.dataretriever.client import QueryResponse
from sunpy.net.tests.strategies import time_attr
from sunpy.time import parse_time


@pytest.fixture
def suvi_client():
    return goes.SUVIClient()


@settings(max_examples=5)
@given(time_attr())
def test_can_handle_query(time):
    suvi_client = goes.SUVIClient()
    assert suvi_client._can_handle_query(time, a.Instrument.suvi) is True
    assert suvi_client._can_handle_query(time, a.Instrument.suvi, a.Wavelength(131 * u.Angstrom)) is True
    assert suvi_client._can_handle_query(time) is False
    assert suvi_client._can_handle_query(time, a.Instrument.aia) is False


def mock_query_object(suvi_client):
    obj = {
        'Start Time': parse_time('2019/05/25 00:50'),
        'End Time': parse_time('2019/05/25 00:52'),
        'Instrument': 'SUVI',
        'Physobs': 'flux',
        'Source': 'GOES',
        'Provider': 'NOAA',
        'Level': '2',
        'Wavelength': 94 * u.Angstrom,
        'url': 'https://mock.url/suvi-l2-ci094.fits'
    }
    return QueryResponse([obj], client=suvi_client)


@pytest.mark.remote_data
@pytest.mark.parametrize(('start', 'end', 'wave', 'level', 'expected_num_files'), [
    ('2019/05/25 00:50', '2019/05/25 00:54', 94, '1b', 6),
    ('2019/05/25 00:50', '2019/05/25 00:54', 304, '2', 1),
])
def test_combined_search(suvi_client, start, end, wave, level, expected_num_files):
    goes_sat = a.goes.SatelliteNumber.sixteen
    qresponse = suvi_client.search(a.Time(start, end), a.Wavelength(wave * u.Angstrom), goes_sat, a.Level(level))
    print("Query Response Length:", len(qresponse))
    # print("Query Response Table:")
    # print(qresponse.table)
    print(f"Start: {start}, End: {end}, Wavelength: {wave}, Level: {level}")
    assert len(qresponse) == expected_num_files



@pytest.mark.remote_data
def test_get_all_wavelengths_level2(suvi_client):
    """Check retrieval for all wavelengths without specifying one."""
    qresponse = suvi_client.search(a.Time('2019/05/25 00:50', '2019/05/25 00:52'),
                                   a.goes.SatelliteNumber.sixteen, a.Level(2))
    assert len(qresponse) == 6


@pytest.mark.remote_data
def test_fetch_working_mock(suvi_client, mocker):
    mocker.patch.object(suvi_client, 'search', return_value=mock_query_object(suvi_client))
    mocker.patch.object(suvi_client, 'fetch', return_value=['mockfile.fits'])

    tr = a.Time('2019/05/25 00:50', '2019/05/25 00:52')
    wave = a.Wavelength(94 * u.Angstrom)
    goes_sat = a.goes.SatelliteNumber.sixteen

    qr = suvi_client.search(tr, a.Instrument.suvi, wave, goes_sat, a.Level(2))
    download_list = suvi_client.fetch(qr)

    assert len(download_list) == len(qr)


def test_show(suvi_client):
    mock_qr = mock_query_object(suvi_client)
    qrshow = mock_qr.show()
    assert {'Start Time', 'End Time', 'Instrument', 'Physobs', 'Source', 'Provider', 'Level', 'Wavelength', 'url'} \
           .issubset(set(qrshow.colnames))
