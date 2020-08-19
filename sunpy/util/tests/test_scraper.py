import datetime
from unittest.mock import Mock, patch

import pytest

import astropy.units as u
from astropy.time import TimeDelta

from sunpy.data.test import rootdir
from sunpy.time import TimeRange, parse_time
from sunpy.util.scraper import Scraper, get_timerange_from_exdict

PATTERN_EXAMPLES = [
    ('%b%y', TimeDelta(31*u.day)),
    ('%m%y', TimeDelta(31*u.day)),
    ('%H%d', TimeDelta(1*u.hour)),
    ('%y%b', TimeDelta(31*u.day)),
]


def testDirectoryDatePattern():
    s = Scraper('%Y/%m/%d/%Y%m%d_%H%M%S_59.fit.gz')
    testpath = '2014/03/05/20140305_013000_59.fit.gz'
    d = parse_time((2014, 3, 5, 1, 30))
    assert s.matches(testpath, d)


def testDirectoryDatePatternFalse():
    s = Scraper('%Y/%m/%d/%Y%m%d_%H%M%S_59.fit.gz')
    testpath = '2013/03/05/20140305_013000_59.fit.gz'
    d = parse_time((2014, 3, 5, 1, 30))
    assert not s.matches(testpath, d)


def testDirectoryObsPattern():
    s = Scraper('%y%m%d/{observatory}_%Y%m%d.fits', observatory='SDO')
    testpath = '140305/SDO_20140305.fits'
    d = parse_time((2014, 3, 5))
    assert s.matches(testpath, d)


def testDirectoryRange():
    s = Scraper('%Y/%m/%d/%Y%m%d_%H.fit.gz')
    directory_list = ['2009/12/30/', '2009/12/31/', '2010/01/01/',
                      '2010/01/02/', '2010/01/03/']
    timerange = TimeRange('2009-12-30', '2010-01-03')
    assert s.range(timerange) == directory_list


def testDirectoryRangeFalse():
    s = Scraper('%Y%m%d/%Y%m%d_%H.fit.gz')
    directory_list = ['20091230/', '20091231/', '20100101/',
                      '20090102/', '20090103/']
    timerange = TimeRange('2009/12/30', '2010/01/03')
    assert s.range(timerange) != directory_list


def testNoDateDirectory():
    s = Scraper('mySpacecraft/myInstrument/xMinutes/aaa%y%b.ext')
    directory_list = ['mySpacecraft/myInstrument/xMinutes/']
    timerange = TimeRange('2009/11/20', '2010/01/03')
    assert s.range(timerange) == directory_list


@pytest.mark.parametrize('pattern, mintime', PATTERN_EXAMPLES)
def test_smallerPattern(pattern, mintime):
    assert mintime == Scraper('')._smallerPattern(pattern)


def testDirectoryRangeHours():
    s = Scraper('%Y%m%d_%H/%H%M.csv')
    timerange = TimeRange('2009-12-31T23:40:00', '2010-01-01T01:15:00')
    assert len(s.range(timerange)) == 3  # 3 directories (1 per hour)


def testDirectoryRange_single():
    s = Scraper('%Y%m%d/%H_%M.csv')
    startdate = parse_time((2010, 10, 10, 5, 0))
    enddate = parse_time((2010, 10, 10, 7, 0))
    timerange = TimeRange(startdate, enddate)
    assert len(s.range(timerange)) == 1


def testDirectoryRange_Month():
    s = Scraper('%Y%m/%d/%j_%H.txt')
    startdate = parse_time((2008, 2, 20, 10))
    enddate = parse_time((2008, 3, 2, 5))
    timerange = TimeRange(startdate, enddate)
    assert len(s.range(timerange)) == 12
    startdate = parse_time((2009, 2, 20, 10))
    enddate = parse_time((2009, 3, 2, 5))
    timerange = TimeRange(startdate, enddate)
    assert len(s.range(timerange)) == 11


def testNoDirectory():
    s = Scraper('files/%Y%m%d_%H%M.dat')
    startdate = parse_time((2010, 1, 10, 20, 30))
    enddate = parse_time((2010, 1, 20, 20, 30))
    timerange = TimeRange(startdate, enddate)
    assert len(s.range(timerange)) == 1


def testExtractDates_usingPattern():
    # Standard pattern
    s = Scraper('data/%Y/%m/%d/fits/swap/swap_00174_fd_%Y%m%d_%H%M%S.fts.gz')
    testURL = 'data/2014/05/14/fits/swap/swap_00174_fd_20140514_200135.fts.gz'
    timeURL = parse_time((2014, 5, 14, 20, 1, 35))
    assert s._extractDateURL(testURL) == timeURL
    # Not-full repeated pattern
    s = Scraper('data/%Y/fits/swap/swap_00174_fd_%Y%m%d_%H%M%S.fts.gz')
    testURL = 'data/2014/fits/swap/swap_00174_fd_20140514_200135.fts.gz'
    timeURL = parse_time((2014, 5, 14, 20, 1, 35))
    assert s._extractDateURL(testURL) == timeURL


def testExtractDates_notSeparators():
    s = Scraper('data/%Y/%m/swap%m%d_%H%M%S')
    testURL = 'data/2014/05/swap0514_200135'
    timeURL = parse_time((2014, 5, 14, 20, 1, 35))
    assert s._extractDateURL(testURL) == timeURL


def testExtractDates_notSeparators_andSimilar():
    s = Scraper('data/%Y/Jun%b%d_%H%M%S')
    testURL = 'data/2014/JunJun14_200135'
    timeURL = parse_time((2014, 6, 14, 20, 1, 35))
    assert s._extractDateURL(testURL) == timeURL
    testURL = 'data/2014/JunMay14_200135'
    timeURL = parse_time((2014, 5, 14, 20, 1, 35))
    assert s._extractDateURL(testURL) == timeURL
    # and testing with the month afterwards
    s = Scraper('data/%Y/%dJun%b_%H%M%S')
    testURL = 'data/2014/14JunJun_200135'
    timeURL = parse_time((2014, 6, 14, 20, 1, 35))
    assert s._extractDateURL(testURL) == timeURL


def testURL_pattern():
    s = Scraper('fd_%Y%m%d_%H%M%S.fts')
    assert s._URL_followsPattern('fd_20130410_231211.fts')
    assert not s._URL_followsPattern('fd_20130410_231211.fts.gz')
    assert not s._URL_followsPattern('fd_20130410_ar_231211.fts.gz')


def testURL_patternMillisecondsGeneric():
    s = Scraper('fd_%Y%m%d_%H%M%S_%e.fts')
    assert s._URL_followsPattern('fd_20130410_231211_119.fts')
    assert not s._URL_followsPattern('fd_20130410_231211.fts.gz')
    assert not s._URL_followsPattern('fd_20130410_ar_231211.fts.gz')


def testURL_patternMillisecondsZeroPadded():
    # Asserts solution to ticket #1954.
    # Milliseconds must be zero-padded in order to match URL lengths.
    now_mock = Mock(return_value=datetime.datetime(2019, 4, 19, 0, 0, 0, 4009))
    with patch('datetime.datetime', now=now_mock):
        s = Scraper('fd_%Y%m%d_%H%M%S_%e.fts')
    now_mock.assert_called_once()
    assert s.now == 'fd_20190419_000000_004.fts'


def testFilesRange_sameDirectory_local():
    s = Scraper('/'.join(['file:/', rootdir,
                          'EIT', 'efz%Y%m%d.%H%M%S_s.fits']))
    startdate = parse_time((2004, 3, 1, 4, 0))
    enddate = parse_time((2004, 3, 1, 6, 30))
    assert len(s.filelist(TimeRange(startdate, enddate))) == 3
    startdate = parse_time((2010, 1, 10, 20, 30))
    enddate = parse_time((2010, 1, 20, 20, 30))
    assert len(s.filelist(TimeRange(startdate, enddate))) == 0


@pytest.mark.remote_data
def testFilesRange_sameDirectory_remote():
    pattern = ('http://solarmonitor.org/data/%Y/%m/%d/'
               'fits/{instrument}/'
               '{instrument}_00174_fd_%Y%m%d_%H%M%S.fts.gz')
    s = Scraper(pattern, instrument='swap')
    startdate = parse_time((2014, 5, 14, 0, 0))
    enddate = parse_time((2014, 5, 14, 6, 30))
    timerange = TimeRange(startdate, enddate)
    assert len(s.filelist(timerange)) == 2
    startdate = parse_time((2014, 5, 14, 21, 0))
    enddate = parse_time((2014, 5, 14, 23, 30))
    timerange = TimeRange(startdate, enddate)
    assert len(s.filelist(timerange)) == 0


@pytest.mark.remote_data
def testFilesRange_sameDirectory_months_remote():
    pattern = ('http://www.srl.caltech.edu/{spacecraft}/DATA/{instrument}/'
               'Ahead/1minute/AeH%y%b.1m')
    s = Scraper(pattern, spacecraft='STEREO', instrument='HET')
    startdate = parse_time((2007, 8, 1))
    enddate = parse_time((2007, 9, 10))
    timerange = TimeRange(startdate, enddate)
    assert len(s.filelist(timerange)) == 2


@pytest.mark.remote_data
def test_ftp():
    pattern = 'ftp://solar-pub.nao.ac.jp/pub/nsro/norh/data/tcx/%Y/%m/tca%y%m%d'
    s = Scraper(pattern)
    timerange = TimeRange('2016/5/18 15:28:00', '2016/5/20 16:30:50')
    urls = s.filelist(timerange)
    assert urls[0] == ('ftp://solar-pub.nao.ac.jp'
                       '/pub/nsro/norh/data/tcx/2016/05/tca160519')
    assert len(urls) == 2


@pytest.mark.remote_data
def test_filelist_url_missing_directory():
    # Asserts solution to ticket #2684.
    # Attempting to access data for the year 1960 results in a 404, so no files are returned.
    pattern = 'http://lasp.colorado.edu/eve/data_access/evewebdataproducts/level2/%Y/%j/'
    s = Scraper(pattern)
    timerange = TimeRange('1960/01/01 00:00:00', '1960/01/02 00:00:00')
    assert len(s.filelist(timerange)) == 0


@pytest.mark.remote_data
def test_filelist_relative_hrefs():
    # the url opened by the scraper from below pattern contains some links which don't have hrefs
    pattern = 'http://www.bbso.njit.edu/pub/archive/%Y/%m/%d/bbso_halph_fr_%Y%m%d_%H%M%S.fts'
    s = Scraper(pattern)
    timerange = TimeRange('2016/5/18 15:28:00', '2016/5/18 16:30:00')
    assert s.domain == 'http://www.bbso.njit.edu/'
    # hrefs are relative to domain here, not to the directory they are present in
    # this checks that `scraper.filelist` returns fileurls relative to the domain
    fileurls = s.filelist(timerange)
    assert fileurls[1] == s.domain + 'pub/archive/2016/05/18/bbso_halph_fr_20160518_160033.fts'


@pytest.mark.parametrize('pattern, check_file', [
    (r'MyFile_%Y_%M_%e\.(\D){2}\.fits', 'MyFile_2020_55_234.aa.fits'),
    (r'(\d){5}_(\d){2}\.fts', '01122_25.fts'),
    (r'_%Y%m%d__%ec(\d){5}_(\d){2}\s.fts', '_20201535__012c12345_33 .fts')])
def test_regex(pattern, check_file):
    s = Scraper(pattern, regex=True)
    assert s._URL_followsPattern(check_file)


@pytest.mark.remote_data
def test_regex_data():
    prefix = r'https://gong2.nso.edu/oQR/zqs/'
    pattern = prefix + r'%Y%m/mrzqs%y%m%d/mrzqs%y%m%dt%H%Mc(\d){4}_(\d){3}\.fits.gz'
    s = Scraper(pattern, regex=True)
    timerange = TimeRange('2020-01-05', '2020-01-06T16:00:00')
    assert s._URL_followsPattern(prefix + '202001/mrzqs200106/mrzqs200106t1514c2226_297.fits.gz')
    assert len(s.filelist(timerange)) == 37


@pytest.mark.remote_data
def test_extract_files_meta():
    baseurl0 = r'ftp://solar-pub.nao.ac.jp/pub/nsro/norh/data/tcx/%Y/%m/(\w){3}%y%m%d'
    extractpattern0 = '{}/tcx/{year:4d}/{month:2d}/{wave}{:4d}{day:2d}'
    s0 = Scraper(baseurl0, regex=True)
    timerange0 = TimeRange('2020/1/1 4:00', '2020/1/2')
    matchdict = {'wave': ['tca', 'tcz']}
    metalist0 = s0._extract_files_meta(timerange0, extractpattern0, matcher=matchdict)
    assert metalist0[0]['wave'] == 'tca'
    assert metalist0[3]['wave'] == 'tcz'
    assert metalist0[1]['day'] == 2

    prefix = r'https://gong2.nso.edu/oQR/zqs/'
    baseurl1 = prefix + r'%Y%m/mrzqs%y%m%d/mrzqs%y%m%dt%H%Mc(\d){4}_(\d){3}\.fits.gz'
    extractpattern1 = ('{}/zqs/{year:4d}{month:2d}/mrzqs{:4d}{day:2d}/mrzqs{:6d}t'
                       '{hour:2d}{minute:2d}c{CAR_ROT:4d}_{:3d}.fits.gz')
    s1 = Scraper(baseurl1, regex=True)
    timerange1 = TimeRange('2020-01-05', '2020-01-05T16:00:00')
    metalist1 = s1._extract_files_meta(timerange1, extractpattern1)
    urls = s1.filelist(timerange1)
    assert metalist1[3]['CAR_ROT'] == 2226
    assert metalist1[-1]['url'] == urls[-1]


@pytest.mark.parametrize('exdict, start, end', [
    ({"year": 2000}, '2000-01-01 00:00:00', '2000-12-31 23:59:59.999000'),
    ({"year": 2016, "month": 2}, '2016-02-01 00:00:00', '2016-02-29 23:59:59.999000'),
    ({'year': 2019, 'month': 2, 'day': 28}, '2019-02-28 00:00:00', '2019-02-28 23:59:59.999000'),
    ({'year': 2020, 'month': 7, 'day': 31, 'hour': 23, 'minute': 59, 'second': 59},
     '2020-07-31 23:59:59', '2020-07-31 23:59:59.999000')])
def test_get_timerange_with_extractor(exdict, start, end):
    tr = TimeRange(start, end)
    file_timerange = get_timerange_from_exdict(exdict)
    assert file_timerange == tr
