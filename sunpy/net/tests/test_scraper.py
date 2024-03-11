import datetime
from pathlib import Path
from unittest.mock import Mock, patch

import pytest

from sunpy.data.test import rootdir
from sunpy.extern import parse
from sunpy.net.scraper import Scraper
from sunpy.time import TimeRange, parse_time


def testDirectoryDatePattern():
    s = Scraper('{{year:4d}}/{{month:2d}}/{{day:2d}}/{{year:4d}}{{month:2d}}{{day:2d}}_{{hour:2d}}{{minute:2d}}{{second:2d}}_59.fit.gz')
    testpath = str(Path('2014/03/05/') / '20140305_013000_59.fit.gz')
    d = parse_time((2014, 3, 5, 1, 30))
    assert s.matches(testpath, d)


def testDirectoryDatePatternFalse():
    s = Scraper('{{year:4d}}/{{month:2d}}/{{day:2d}}/{{year:4d}}{{month:2d}}{{day:2d}}_{{hour:2d}}{{minute:2d}}{{second:2d}}_59.fit.gz')
    testpath = str(Path('2013/03/05/') / '20140305_013000_59.fit.gz')
    d = parse_time((2014, 3, 5, 1, 30))
    assert not s.matches(testpath, d)


def testDirectoryObsPattern():
    s = Scraper('{{year:2d}}{{month:2d}}{{day:2d}}/{observatory}_{{year:4d}}{{month:2d}}{{day:2d}}.fits', observatory='SDO')
    testpath = str(Path('140305') / 'SDO_20140305.fits')
    d = parse_time((2014, 3, 5))
    assert s.matches(testpath, d)


def testDirectoryRange():
    s = Scraper('{{year:4d}}/{{month:2d}}/{{day:2d}}/{{year:4d}}{{month:2d}}{{day:2d}}_{{hour:2d}}.fit.gz')
    directory_list = ['2009/12/30/', '2009/12/31/', '2010/01/01/',
                      '2010/01/02/', '2010/01/03/']
    timerange = TimeRange('2009-12-30', '2010-01-03')
    assert s.range(timerange) == directory_list


def testDirectoryRegex():
    # Test for Windows where '\' is a path separator and not part of the regex
    s = Scraper('scheme://a.url.with/a/few/forward/slashes/andbacklash\\inthename.ext', regex=True)
    timerange = TimeRange('2019-02-01', '2019-02-03')
    directory = s.range(timerange)
    assert directory == ['scheme://a.url.with/a/few/forward/slashes/']


def testDirectoryRangeFalse():
    s = Scraper('{{year:4d}}{{month:2d}}{{day:2d}}/{{year:4d}}{{month:2d}}{{day:2d}}_{{hour:2d}}.fit.gz')
    directory_list = ['20091230/', '20091231/', '20100101/',
                      '20090102/', '20090103/']
    timerange = TimeRange('2009/12/30', '2010/01/03')
    assert s.range(timerange) != directory_list


def testNoDateDirectory():
    s = Scraper('mySpacecraft/myInstrument/xMinutes/aaa{{year:4d}}{{month_name_abbr:l}}.ext')
    directory_list = ['mySpacecraft/myInstrument/xMinutes/']
    timerange = TimeRange('2009/11/20', '2010/01/03')
    assert s.range(timerange) == directory_list

def testDirectoryRangeHours():
    s = Scraper('{{year:4d}}{{month:2d}}{{day:2d}}_{{hour:2d}}/{{hour:2d}}{{minute:2d}}.csv')
    timerange = TimeRange('2009-12-31T23:40:00', '2010-01-01T01:15:00')
    assert len(s.range(timerange)) == 3  # 3 directories (1 per hour)


def testDirectoryRange_single():
    s = Scraper('{{year:4d}}{{month:2d}}{{day:2d}}/{{hour:2d}}_{{minute:2d}}.csv')
    startdate = parse_time((2010, 10, 10, 5, 0))
    enddate = parse_time((2010, 10, 10, 7, 0))
    timerange = TimeRange(startdate, enddate)
    assert len(s.range(timerange)) == 1


def testDirectoryRange_Month():
    s = Scraper('{{year:4d}}{{month:2d}}/{{day:2d}}/{{day_of_year:3d}}_{{hour:2d}}.txt')
    startdate = parse_time((2008, 2, 20, 10))
    enddate = parse_time((2008, 3, 2, 5))
    timerange = TimeRange(startdate, enddate)
    assert len(s.range(timerange)) == 12
    startdate = parse_time((2009, 2, 20, 10))
    enddate = parse_time((2009, 3, 2, 5))
    timerange = TimeRange(startdate, enddate)
    assert len(s.range(timerange)) == 11


def testNoDirectory():
    s = Scraper('files/{{year:4d}}{{month:2d}}{{day:2d}}_{{hour:2d}}{{month:2d}}.dat')
    startdate = parse_time((2010, 1, 10, 20, 30))
    enddate = parse_time((2010, 1, 20, 20, 30))
    timerange = TimeRange(startdate, enddate)
    assert len(s.range(timerange)) == 1

@pytest.mark.parametrize(('pattern', 'filename', 'metadict'), [
    ('_{{year:4d}}{{month:2d}}{{day:2d}}__{{millisecond:3d}}c{{:5d}}_{{:2d}}{{}}.fts', '_20201535__012c12345_33 .fts',
     {'year': 2020, 'month': 15, 'day': 35, 'millisecond': 12}),
    ('fd_{{year:4d}}{{month:2d}}{{day:2d}}_{{hour:2d}}{{minute:2d}}{{second:2d}}_{{millisecond:3d}}.fts', 'fd_20130410_231211_119.fts',
     {'year': 2013, 'month': 4, 'day': 10, 'hour': 23, 'minute': 12, 'second': 11, 'millisecond': 119})
])
def testURL_pattern(pattern, filename, metadict):
    assert parse(pattern.format(None), filename).named == metadict

def testURL_patternMillisecondsZeroPadded():
    # Asserts solution to ticket #1954.
    # Milliseconds must be zero-padded in order to match URL lengths.
    now_mock = Mock(return_value=datetime.datetime(2019, 4, 19, 0, 0, 0, 4009))
    with patch('sunpy.net.scraper.datetime', now=now_mock):
        s = Scraper('fd_{{year:4d}}{{month:2d}}{{day:2d}}_{{hour:2d}}{{minute:2d}}{{second:2d}}_{{millisecond:3d}}.fts')
    now_mock.assert_called_once()
    assert s.now == 'fd_20190419_000000_004.fts'


def testFilesRange_sameDirectory_local():
    s = Scraper('/'.join(['file:/', str(rootdir),
                          'EIT_header', 'efz{{year:4d}}{{month:2d}}{{day:2d}}.{{hour:2d}}{{minute:2d}}{{second:2d}}_s.header']))
    startdate = parse_time((2004, 3, 1, 4, 0))
    enddate = parse_time((2004, 3, 1, 6, 30))
    assert len(s.filelist(TimeRange(startdate, enddate))) == 3
    startdate = parse_time((2010, 1, 10, 20, 30))
    enddate = parse_time((2010, 1, 20, 20, 30))
    assert len(s.filelist(TimeRange(startdate, enddate))) == 0


@pytest.mark.remote_data
def testFilesRange_sameDirectory_remote():
    pattern = ('http://proba2.oma.be/{instrument}/data/bsd/{{year:4d}}/{{month:2d}}/{{day:2d}}/'
               '{instrument}_lv1_{{year:4d}}{{month:2d}}{{day:2d}}_{{hour:2d}}{{minute:2d}}{{second:2d}}.fits')
    s = Scraper(pattern, instrument='swap')
    startdate = parse_time((2014, 5, 14, 0, 0))
    enddate = parse_time((2014, 5, 14, 0, 5))
    timerange = TimeRange(startdate, enddate)
    assert len(s.filelist(timerange)) == 2
    startdate = parse_time((2014, 5, 14, 0, 6))
    enddate = parse_time((2014, 5, 14, 0, 7))
    timerange = TimeRange(startdate, enddate)
    assert len(s.filelist(timerange)) == 0


@pytest.mark.remote_data
def testFilesRange_sameDirectory_months_remote():
    pattern = ('http://www.srl.caltech.edu/{spacecraft}/DATA/{instrument}/'
               'Ahead/1minute/AeH{{year:2d}}{{month_name_abbr:w}}.1m')
    s = Scraper(pattern, spacecraft='STEREO', instrument='HET')
    startdate = parse_time((2007, 8, 1))
    enddate = parse_time((2007, 9, 10))
    timerange = TimeRange(startdate, enddate)
    files = s.filelist(timerange)
    assert files == ['http://www.srl.caltech.edu/STEREO/DATA/HET/Ahead/1minute/AeH07Aug.1m',
                     'http://www.srl.caltech.edu/STEREO/DATA/HET/Ahead/1minute/AeH07Sep.1m']


@pytest.mark.remote_data
def test_ftp():
    pattern = 'ftp://solar-pub.nao.ac.jp/pub/nsro/norh/data/tcx/{{year:4d}}/{{month:2d}}/tca{{year:4d}}{{month:2d}}{{day:2d}}'
    s = Scraper(pattern)
    timerange = TimeRange('2016/5/18 15:28:00', '2016/5/20 16:30:50')
    urls = s.filelist(timerange)
    assert urls[0] == ('ftp://solar-pub.nao.ac.jp'
                       '/pub/nsro/norh/data/tcx/2016/05/tca160518')
    assert len(urls) == 3


@pytest.mark.remote_data
def test_filelist_url_missing_directory():
    # Asserts solution to ticket #2684.
    # Attempting to access data for the year 1960 results in a 404, so no files are returned.
    pattern = 'http://lasp.colorado.edu/eve/data_access/evewebdataproducts/level2/{{year:4d}}/{{day_of_year:3d}}/'
    s = Scraper(pattern)
    timerange = TimeRange('1960/01/01 00:00:00', '1960/01/02 00:00:00')
    assert len(s.filelist(timerange)) == 0


@pytest.mark.remote_data
def test_filelist_relative_hrefs():
    # the url opened by the scraper from below pattern contains some links which don't have hrefs
    pattern = 'http://www.bbso.njit.edu/pub/archive/{{year:4d}}/{{month:2d}}/{{day:2d}}/bbso_halph_fr_{{year:4d}}{{month:2d}}{{day:2d}}_{{hour:2d}}{{minute:2d}}{{second:2d}}.fts'
    s = Scraper(pattern)
    timerange = TimeRange('2016/5/18 15:28:00', '2016/5/18 16:30:00')
    assert s.domain == 'http://www.bbso.njit.edu/'
    # hrefs are relative to domain here, not to the directory they are present in
    # this checks that `scraper.filelist` returns fileurls relative to the domain
    fileurls = s.filelist(timerange)
    assert fileurls[1] == s.domain + 'pub/archive/2016/05/18/bbso_halph_fr_20160518_160033.fts'


@pytest.mark.remote_data
def test_parse_pattern_data():
    prefix = 'https://gong2.nso.edu/oQR/zqs/'
    pattern = prefix + '{{year:4d}}{{month:2d}}/mrzqs{{year:2d}}{{month:2d}}{{day:2d}}/mrzqs{{year:2d}}{{month:2d}}{{day:2d}}t{{hour:2d}}{{minute:2d}}c{{:4d}}_{{:3d}}.fits.gz'
    s = Scraper(pattern)
    timerange = TimeRange('2020-01-05', '2020-01-06T16:00:00')
    assert len(s.filelist(timerange)) == 37


@pytest.mark.remote_data
def test_extract_files_meta():
    pattern0 = 'ftp://solar-pub.nao.ac.jp/pub/nsro/norh/data/tcx/{{year:4d}}/{{month:2d}}/{{wave}}{{year:2d}}{{month:2d}}{{day:2d}}'
    s0 = Scraper(pattern0)
    timerange0 = TimeRange('2020/1/1 4:00', '2020/1/2')
    matchdict = {'wave': ['tca', 'tcz']}
    metalist0 = s0._extract_files_meta(timerange0, matcher=matchdict)
    assert metalist0[0]['wave'] == 'tca'
    assert metalist0[3]['wave'] == 'tcz'
    assert metalist0[1]['day'] == 2

    pattern1 = ('https://gong2.nso.edu/oQR/zqs/{{year:4d}}{{month:2d}}/mrzqs{{year:2d}}{{month:2d}}{{day:2d}}'
                '/mrzqs{{year:2d}}{{month:2d}}{{day:2d}}t{{hour:2d}}{{minute:2d}}c{{CAR_ROT:4d}}_{{:3d}}.fits.gz')
    s1 = Scraper(pattern1)
    timerange1 = TimeRange('2020-01-05', '2020-01-05T16:00:00')
    metalist1 = s1._extract_files_meta(timerange1)
    urls = s1.filelist(timerange1)
    assert metalist1[3]['CAR_ROT'] == 2226
    assert metalist1[-1]['url'] == urls[-1]

@pytest.mark.remote_data
def test_yearly_overlap():
    # Check that a time range that falls within the interval that a file represents
    # returns a single result.
    pattern = "https://www.ngdc.noaa.gov/stp/space-weather/solar-data/solar-features/solar-flares/x-rays/goes/xrs/goes-xrs-report_{{year:4d}}.txt"
    scraper = Scraper(pattern)

    # Should return a single file for 2013
    trange = TimeRange("2013-01-02", "2013-01-03")
    assert len(scraper.filelist(trange)) == 1
