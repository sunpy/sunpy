import tempfile

import pytest
from hypothesis import given

import astropy.units as u
from astropy.time import Time

import sunpy.net.dataretriever.sources.goes as goes
from sunpy.net import Fido
from sunpy.net import attrs as a
from sunpy.net.dataretriever.client import QueryResponse
from sunpy.net.tests.strategies import time_attr
from sunpy.time import parse_time


@pytest.fixture
def suvi_client():
    return goes.SUVIClient()


@given(time_attr())
def test_can_handle_query(time):
    # Don't use the fixture, as hypothesis complains
    suvi_client = goes.SUVIClient()
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


def mock_query_object(suvi_client):
    """
    Creating a Query Response object and prefilling it with some information
    """
    # Creating a Query Response Object
    start = '2019/05/25 00:50'
    end = '2019/05/25 00:52'
    wave = 94 * u.Angstrom
    obj = {
        'Start Time': parse_time(start),
        'End Time': parse_time(end),
        'Instrument': 'SUVI',
        'Physobs': 'flux',
        'Source': 'GOES',
        'Provider': 'NOAA',
        'Level': '2',
        'Wavelength': wave,
        'url': ('https://data.ngdc.noaa.gov/platforms/solar-space-observing-satellites'
                '/goes/goes16/l2/data/suvi-l2-ci094/2019/05/25/'
                'dr_suvi-l2-ci094_g16_s20190525T005200Z_e20190525T005600Z_v1-0-0.fits')
    }
    results = QueryResponse([obj], client=suvi_client)
    return results


def test_attr_reg():
    a.Instrument.suvi = a.Instrument("SUVI")
    a.goes.SatelliteNumber.A16 = a.goes.SatelliteNumber("16")


@pytest.mark.remote_data
def test_fetch_working(suvi_client):
    """
    Tests if the online server for goes_suvi is working.
    This also checks if the mock is working well.
    """
    start = '2019/05/25 00:50'
    end = '2019/05/25 00:52'
    wave = 94 * u.Angstrom
    goes_sat = a.goes.SatelliteNumber.sixteen
    tr = a.Time(start, end)
    qr1 = suvi_client.search(tr, a.Instrument.suvi, a.Wavelength(wave), goes_sat, a.Level(2))

    # Mock QueryResponse object
    mock_qr = mock_query_object(suvi_client)

    # Compare if two objects have the same attribute

    mock_qr = mock_qr[0]
    qr = qr1[0]

    assert mock_qr['Source'] == qr['Source']
    assert mock_qr['Provider'] == qr['Provider']
    assert mock_qr['Physobs'] == qr['Physobs']
    assert mock_qr['Instrument'] == qr['Instrument']
    assert mock_qr['url'] == qr['url']

    assert qr1['Start Time'] == Time("2019-05-25T00:52:00.000")
    assert qr1['End Time'] == Time("2019-05-25T00:56:00.000")

    with tempfile.TemporaryDirectory() as tmpdirname:
        download_list = suvi_client.fetch(qr1, path=tmpdirname)
    assert len(download_list) == len(qr1)


@pytest.mark.remote_data
@pytest.mark.parametrize(('start', 'end', 'wave', 'expected_num_files'),
                         [('2019/05/25 00:50', '2019/05/25 00:52', 94, 1),
                          ('2019/05/25 00:50', '2019/05/25 00:52', 131, 1),
                          ('2019/05/25 00:50', '2019/05/25 00:52', 171, 1),
                          ('2019/05/25 00:50', '2019/05/25 00:52', 195, 1),
                          ('2019/05/25 00:50', '2019/05/25 00:52', 284, 1),
                          ('2019/05/25 00:50', '2019/05/25 00:52', 304, 1)]
                         )
def test_get_url_for_time_range_level2(suvi_client, start, end, wave, expected_num_files):
    goes_sat = a.goes.SatelliteNumber.sixteen
    qresponse = suvi_client.search(a.Time(start, end), a.Wavelength(wave * u.Angstrom), goes_sat, a.Level(2))
    urls = [i['url'] for i in qresponse]
    assert isinstance(urls, list)
    assert len(urls) == expected_num_files


@pytest.mark.remote_data
@pytest.mark.parametrize(('start', 'end', 'expected_num_files'),
                         [('2019/05/25 00:50', '2019/05/25 00:52', 6)]
                         )
def test_get_url_for_time_range_level2_allwave(suvi_client, start, end, expected_num_files):
    """check that we get all wavelengths if no wavelength is given"""
    goes_sat = a.goes.SatelliteNumber.sixteen
    qresponse = suvi_client.search(a.Time(start, end), goes_sat, a.Level(2))
    urls = [i['url'] for i in qresponse]
    assert isinstance(urls, list)
    assert len(urls) == expected_num_files


@pytest.mark.remote_data
@pytest.mark.parametrize(('start', 'end', 'wave', 'expected_num_files'),
                         [('2019/05/25 00:50', '2019/05/25 00:54', 94, 6),
                          ('2019/05/25 00:50', '2019/05/25 00:54', 131, 3),
                          ('2019/05/25 00:50', '2019/05/25 00:54', 171, 2),
                          ('2019/05/25 00:50', '2019/05/25 00:54', 195, 7),
                          ('2019/05/25 00:50', '2019/05/25 00:54', 284, 2),
                          ('2019/05/25 00:50', '2019/05/25 00:54', 304, 4)]
                         )
def test_get_url_for_time_range_level1b(suvi_client, start, end, wave, expected_num_files):
    """check that we get all wavelengths if no wavelength is given"""
    goes_sat = a.goes.SatelliteNumber.sixteen
    qresponse = suvi_client.search(a.Time(start, end), a.Wavelength(
        wave * u.Angstrom), goes_sat, a.Level('1b'))
    urls = [i['url'] for i in qresponse]
    assert isinstance(urls, list)
    assert len(urls) == expected_num_files


@pytest.mark.remote_data
@pytest.mark.parametrize(('start', 'end', 'wave', 'expected_num_files'),
                         [('2019/05/25 00:50', '2019/05/25 00:54', 94, 6),
                          ('2019/05/25 00:50', '2019/05/25 00:54', 131, 3),
                          ('2019/05/25 00:50', '2019/05/25 00:54', 171, 2),
                          ('2019/05/25 00:50', '2019/05/25 00:54', 195, 7),
                          ('2019/05/25 00:50', '2019/05/25 00:54', 284, 2),
                          ('2019/05/25 00:50', '2019/05/25 00:54', 304, 4)]
                         )
def test_fido_onewave_level1b(start, end, wave, expected_num_files):
    goes_sat = a.goes.SatelliteNumber.sixteen
    result = Fido.search(a.Time(start, end), a.Instrument.suvi, goes_sat,
                         a.Wavelength(wave * u.Angstrom), a.Level('1b'))
    assert result.file_num == expected_num_files


@pytest.mark.remote_data
@pytest.mark.parametrize(('start', 'end', 'wave1', 'wave2', 'expected_num_files'),
                         [('2019/05/25 00:50', '2019/05/25 00:54', 1, 100, 6),
                          ('2019/05/25 00:50', '2019/05/25 00:54', 1, 150, 9),
                          ('2019/05/25 00:50', '2019/05/25 00:54', 1, 180, 11),
                          ('2019/05/25 00:50', '2019/05/25 00:54', 1, 200, 18),
                          ('2019/05/25 00:50', '2019/05/25 00:54', 1, 300, 20),
                          ('2019/05/25 00:50', '2019/05/25 00:54', 1, 310, 24)]
                         )
def test_fido_waverange_level1b(start, end, wave1, wave2, expected_num_files):
    """check that we get all wavelengths if no wavelength is given"""
    goes_sat = a.goes.SatelliteNumber.sixteen
    result = Fido.search(a.Time(start, end), a.Instrument.suvi, goes_sat,
                         a.Wavelength(wave1 * u.Angstrom, wave2 * u.Angstrom),
                         a.Level('1b'))
    assert result.file_num == expected_num_files


@pytest.mark.remote_data
@pytest.mark.parametrize(('start', 'end', 'expected_num_files'),
                         [('2019/05/25 00:50', '2019/05/25 00:52', 6)]
                         )
def test_query(suvi_client, start, end, expected_num_files):
    goes_sat = a.goes.SatelliteNumber.sixteen
    qr1 = suvi_client.search(a.Time(start, end), a.Instrument.suvi, goes_sat, a.Level.two)
    assert isinstance(qr1, QueryResponse)
    assert len(qr1) == expected_num_files
    assert qr1['Start Time'][0] == parse_time('2019/05/25 00:52')
    assert qr1['End Time'][1] == parse_time('2019/05/25 00:56')


def test_show(suvi_client):
    mock_qr = mock_query_object(suvi_client)
    qrshow0 = mock_qr.show()
    qrshow1 = mock_qr.show('Start Time', 'Instrument')
    allcols = {'Start Time', 'End Time', 'Instrument', 'Physobs', 'Source',
               'Provider', 'Level', 'Wavelength', 'url'}
    assert not allcols.difference(qrshow0.colnames)
    assert qrshow1.colnames == ['Start Time', 'Instrument']
    assert qrshow0['Instrument'][0] == 'SUVI'
