import io
from unittest import mock
from urllib.error import HTTPError, URLError

import pytest

import sunpy.net.dataretriever.sources.hinode as hinode
from sunpy.net import attrs as a

DARTS = 'https://data.darts.isas.jaxa.jp/pub/hinode/sot/'
LISTINGS = {
    'level0/2013/07/13/FG/H0000/': ['FG20130713_001332.7.fits'],
    'level1/2011/12/13/FG/H0600/': ['FG20111213_060430.3.fits'],
    'level0/2013/07/13/SP4D/H1000/': ['SP4D20130713_100250.1.fits'],
    'level1hao/2011/12/13/SP3D/': ['20111213_025004/'],
    'level1hao/2011/12/13/SP3D/20111213_025004/': [
        '20111213_025004_stksimg.save', 'SP3D20111213_025004.4C.fits', 'SP3D20111213_025009.5C.fits', 'sp_prep.log',
    ],
    'level2hao/2011/12/13/SP3D/': ['20111213_025004/', '20111213_061504/'],
    'level2hao/2011/12/13/SP3D/20111213_025004/': [
        '20111213_025004.fits', '20111213_025004.sav', 'SP3D20111213_025004.4C.fits.xml', 'thumbnails/',
    ],
    'level2hao/2011/12/13/SP3D/20111213_061504/': ['20111213_061504.fits'],
    'level2.1hao/2011/12/13/SP3D/': ['20111213_025004/'],
    'level2.1hao/2011/12/13/SP3D/20111213_025004/': ['20111213_025004_L2.1.fits', '20111213_025004_disambig_log.txt'],
}
SEARCHES = [
    (hinode.HinodeSOTFGClient, 'FG', '0', '2013-07-13 00:13:32', '2013-07-13 00:13:33',
     ['level0/2013/07/13/FG/H0000/FG20130713_001332.7.fits']),
    (hinode.HinodeSOTFGClient, 'FG', '1', '2011-12-13 06:04:30', '2011-12-13 06:04:32',
     ['level1/2011/12/13/FG/H0600/FG20111213_060430.3.fits']),
    (hinode.HinodeSOTSPClient, 'SP', '0', '2013-07-13 10:02:50', '2013-07-13 10:02:51',
     ['level0/2013/07/13/SP4D/H1000/SP4D20130713_100250.1.fits']),
    (hinode.HinodeSOTSPClient, 'SP', '1', '2011-12-13 02:00', '2011-12-13 03:00',
     ['level1hao/2011/12/13/SP3D/20111213_025004/SP3D20111213_025004.4C.fits',
      'level1hao/2011/12/13/SP3D/20111213_025004/SP3D20111213_025009.5C.fits']),
    (hinode.HinodeSOTSPClient, 'SP', '2', '2011-12-13', '2011-12-14',
     ['level2hao/2011/12/13/SP3D/20111213_025004/20111213_025004.fits',
      'level2hao/2011/12/13/SP3D/20111213_061504/20111213_061504.fits']),
    (hinode.HinodeSOTSPClient, 'SP', '2.1', '2011-12-13', '2011-12-14',
     ['level2.1hao/2011/12/13/SP3D/20111213_025004/20111213_025004_L2.1.fits']),
]


def query(detector, level, start='2011-12-13', end='2011-12-14'):
    return a.Time(start, end), a.Instrument('SOT'), a.Level(level), a.Detector(detector)


@pytest.fixture
def urlopen():
    def fake_urlopen(url):
        if (directory := url.removeprefix(DARTS)) not in LISTINGS:
            raise HTTPError(url, 404, 'Not Found', {}, None)
        return io.BytesIO(''.join(f'<a href="{link}"></a>' for link in LISTINGS[directory]).encode())

    with mock.patch('sunpy.net.scraper.urlopen', side_effect=fake_urlopen) as mocked:
        yield mocked


@pytest.mark.parametrize(
    ('client', 'detector', 'level', 'can_handle'),
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
def test_can_handle_query(client, detector, level, can_handle):
    assert client._can_handle_query(*query(detector, level)) is can_handle


@pytest.mark.parametrize(('client', 'detector', 'level', 'start', 'end', 'files'), SEARCHES)
def test_search(urlopen, client, detector, level, start, end, files):
    qr = client().search(*query(detector, level, start, end))
    assert list(qr['url']) == [DARTS + file for file in files]


def test_scan_requests(urlopen):
    hinode.HinodeSOTSPClient().search(*query('SP', '2'))
    # Two day listings and two scans, not one request for each of the 86401 seconds.
    assert urlopen.call_count == 4


@pytest.mark.parametrize(('client', 'detector', 'level'), [
    (hinode.HinodeSOTFGClient, 'FG', '0'),
    (hinode.HinodeSOTSPClient, 'SP', '2'),
])
def test_search_error(client, detector, level):
    # Fido reports what a search raises in ``UnifiedResponse.errors``.
    with mock.patch('sunpy.net.scraper.urlopen', side_effect=URLError('down')), pytest.raises(URLError):
        client().search(*query(detector, level))


def test_show(urlopen):
    qr = hinode.HinodeSOTSPClient().search(*query('SP', '2'))
    allcols = {'Start Time', 'End Time', 'Instrument', 'Source', 'Provider', 'Level', 'Detector', 'url'}
    assert set(qr.show().colnames) == allcols
    assert qr.show('Start Time', 'Detector').colnames == ['Start Time', 'Detector']
    assert tuple(qr['Instrument', 'Level', 'Detector'][0]) == ('SOT', '2', 'SP')


@pytest.mark.remote_data
@pytest.mark.parametrize(('client', 'detector', 'level', 'start', 'end', 'files'), SEARCHES)
def test_real_archive_urls(client, detector, level, start, end, files):
    qr = client().search(*query(detector, level, start, end))
    assert {DARTS + file for file in files} <= set(qr['url'])
