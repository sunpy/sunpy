from unittest import mock

import pytest

import sunpy.net.dataretriever.sources.hinode as hinode
from sunpy.net import attrs as a
from sunpy.net.dataretriever.client import QueryResponse
from sunpy.time import parse_time


@pytest.fixture
def fg_client():
    return hinode.HinodeSOTFGClient()


@pytest.fixture
def sp_client():
    return hinode.HinodeSOTSPClient()


@pytest.mark.parametrize(
    ('client_cls', 'detector', 'level', 'can_handle'),
    [
        (hinode.HinodeSOTFGClient, 'FG', '0', True),
        (hinode.HinodeSOTFGClient, 'FG', '1', True),
        (hinode.HinodeSOTFGClient, 'SP', '1', False),
        (hinode.HinodeSOTFGClient, 'FG', '2', False),
        (hinode.HinodeSOTSPClient, 'SP', '0', True),
        (hinode.HinodeSOTSPClient, 'SP', '1', True),
        (hinode.HinodeSOTSPClient, 'SP', '2', True),
        (hinode.HinodeSOTSPClient, 'SP', '2.1', True),
        (hinode.HinodeSOTSPClient, 'FG', '1', False),
    ],
)
def test_can_handle_query(client_cls, detector, level, can_handle):
    query = (
        a.Time('2011-12-13 06:00', '2011-12-13 07:00'),
        a.Instrument('SOT'),
        a.Level(level),
        a.hinode.SOTDetector(detector),
    )
    assert client_cls._can_handle_query(*query) is can_handle


@pytest.mark.parametrize(
    ('client_cls', 'detector', 'provider', 'can_handle'),
    [
        (hinode.HinodeSOTFGClient, 'FG', 'DARTS', True),
        (hinode.HinodeSOTFGClient, 'FG', 'LMSAL', True),
        (hinode.HinodeSOTSPClient, 'SP', 'DARTS', True),
        (hinode.HinodeSOTSPClient, 'SP', 'LMSAL', True),
        (hinode.HinodeSOTSPClient, 'SP', 'Fake', False),
    ],
)
def test_can_handle_query_provider(client_cls, detector, provider, can_handle):
    query = (
        a.Time('2011-12-13 06:00', '2011-12-13 07:00'),
        a.Instrument('SOT'),
        a.Level(1),
        a.hinode.SOTDetector(detector),
        a.Provider(provider),
    )
    assert client_cls._can_handle_query(*query) is can_handle


REAL_URL_CASES = [
    (
        'FG',
        '0',
        '2013-07-13 00:13:32',
        '2013-07-13 00:13:33',
        'https://data.darts.isas.jaxa.jp/pub/hinode/sot/level0/2013/07/13/FG/H0000/FG20130713_001332.7.fits',
        'DARTS',
    ),
    (
        'FG',
        '1',
        '2011-12-13 06:04:30',
        '2011-12-13 06:04:32',
        'https://data.darts.isas.jaxa.jp/pub/hinode/sot/level1/2011/12/13/FG/H0600/FG20111213_060430.3.fits',
        'DARTS',
    ),
    (
        'SP',
        '0',
        '2013-07-13 10:02:50',
        '2013-07-13 10:02:51',
        'https://data.darts.isas.jaxa.jp/pub/hinode/sot/level0/2013/07/13/SP4D/H1000/SP4D20130713_100250.1.fits',
        'DARTS',
    ),
    (
        'SP',
        '1',
        '2011-12-13 02:50:04',
        '2011-12-13 02:50:05',
        'https://data.darts.isas.jaxa.jp/pub/hinode/sot/level1hao/2011/12/13/SP3D/20111213_025004/SP3D20111213_025004.4C.fits',
        'DARTS',
    ),
    (
        'SP',
        '2',
        '2011-12-13 02:50:04',
        '2011-12-13 02:50:05',
        'https://data.darts.isas.jaxa.jp/pub/hinode/sot/level2hao/2011/12/13/SP3D/20111213_025004/20111213_025004.fits',
        'DARTS',
    ),
    (
        'SP',
        '2.1',
        '2011-12-13 02:50:04',
        '2011-12-13 02:50:05',
        'https://data.darts.isas.jaxa.jp/pub/hinode/sot/level2.1hao/2011/12/13/SP3D/20111213_025004/20111213_025004_L2.1.fits',
        'DARTS',
    ),
    (
        'SP',
        '2',
        '2011-12-13 02:50:04',
        '2011-12-13 02:50:05',
        'https://sot.lmsal.com/data/sot/level2hao/2011/12/13/SP3D/20111213_025004/20111213_025004.fits',
        'LMSAL',
    ),
]


@pytest.mark.remote_data
@pytest.mark.parametrize(
    ('detector', 'level', 'start', 'end', 'expected_url', 'provider'),
    REAL_URL_CASES,
)
def test_real_archive_urls(fg_client, sp_client, detector, level, start, end, expected_url, provider):
    hinode_client = fg_client if detector == 'FG' else sp_client
    qr = hinode_client.search(
        a.Time(start, end),
        a.Instrument('SOT'),
        a.Level(level),
        a.hinode.SOTDetector(detector),
        a.Provider(provider),
    )

    assert len(qr) >= 1
    assert expected_url in list(qr['url'])


@pytest.mark.parametrize(
    ('client_cls', 'detector', 'level', 'provider'),
    [
        (
            hinode.HinodeSOTFGClient,
            'FG',
            1,
            'LMSAL',
        ),
        (
            hinode.HinodeSOTFGClient,
            'FG',
            1,
            'DARTS',
        ),
        (
            hinode.HinodeSOTSPClient,
            'SP',
            2,
            'LMSAL',
        ),
        (
            hinode.HinodeSOTSPClient,
            'SP',
            2,
            'DARTS',
        ),
    ],
)
def test_search_provider_url(client_cls, detector, level, provider):
    lmsal_url = 'https://sot.lmsal.com/data/sot/'
    darts_url = 'https://data.darts.isas.jaxa.jp/pub/hinode/sot/'
    fake_url = (lmsal_url if provider == 'LMSAL' else darts_url) + 'level' + str(level) + 'hao' + '/somefile.fits'
    filesmeta = [{
        'year': 2011, 'month': 12, 'day': 13,
        'hour': 6, 'minute': 0, 'second': 0,
        'url': fake_url,
    }]
    client = client_cls()
    with mock.patch('sunpy.net.dataretriever.client.Scraper._extract_files_meta', return_value=filesmeta):
        qr = client.search(
            a.Time('2011-12-13 06:00', '2011-12-13 06:15'),
            a.Instrument('SOT'),
            a.Level(level),
            a.hinode.SOTDetector(detector),
            a.Provider(provider),
        )
    assert len(qr) == 1
    assert qr[0]['Provider'] == provider
    assert qr[0]['url'].startswith(lmsal_url if provider == 'LMSAL' else darts_url)


def test_search_unsupported_fg_level(fg_client):
    qr = fg_client.search(
        a.Time('2011-12-13 06:00', '2011-12-13 07:00'),
        a.Instrument('SOT'),
        a.Level(2),
        a.hinode.SOTDetector('FG'),
    )
    assert len(qr) == 0


def test_attr_reg():
    assert hasattr(a, 'hinode')
    assert a.hinode.SOTDetector.fg == a.hinode.SOTDetector('FG')
    assert a.hinode.SOTDetector.sp == a.hinode.SOTDetector('SP')


def mock_query_object(client, detector='FG', level='1', url=None):
    if url is None:
        if detector == 'FG':
            url = ('https://data.darts.isas.jaxa.jp/pub/hinode/sot/level1/'
                   '2011/12/13/FG/H0600/FG20111213_060430.3.fits')
        else:
            url = ('https://data.darts.isas.jaxa.jp/pub/hinode/sot/level1hao/'
                   '2011/12/13/SP3D/20111213_025004/SP3D20111213_025004.4C.fits')

    obj = {
        'Start Time': parse_time('2011-12-13T06:00:00'),
        'End Time': parse_time('2011-12-13T06:00:00'),
        'Instrument': 'SOT',
        'SOTDetector': detector,
        'Level': level,
        'Source': 'Hinode',
        'Provider': 'DARTS',
        'url': url,
    }
    return QueryResponse([obj], client=client)


def test_show(fg_client):
    mock_qr = mock_query_object(fg_client)
    qrshow0 = mock_qr.show()
    qrshow1 = mock_qr.show('Start Time', 'SOTDetector')
    allcols = {'Start Time', 'End Time', 'Instrument', 'SOTDetector', 'Level', 'Source', 'Provider', 'url'}
    assert not allcols.difference(qrshow0.colnames)
    assert qrshow1.colnames == ['Start Time', 'SOTDetector']
    assert qrshow0['Instrument'][0] == 'SOT'
