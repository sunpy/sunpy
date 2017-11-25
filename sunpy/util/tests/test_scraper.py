from __future__ import absolute_import, division, print_function

import pytest
import datetime

from sunpy.time import TimeRange
from sunpy.util.scraper import Scraper

PATTERN_EXAMPLES = [
    ('%b%y', datetime.timedelta(days=31)),
    ('%m%y', datetime.timedelta(days=31)),
    ('%H%d', datetime.timedelta(hours=1)),
    ('%y%b', datetime.timedelta(days=31)),
]


def testDirectoryDatePattern():
    s = Scraper('%Y/%m/%d/%Y%m%d_%H%M%S_59.fit.gz')
    testpath = '2014/03/05/20140305_013000_59.fit.gz'
    d = datetime.datetime(2014, 3, 5, 1, 30)
    assert s.matches(testpath, d)


def testDirectoryDatePatternFalse():
    s = Scraper('%Y/%m/%d/%Y%m%d_%H%M%S_59.fit.gz')
    testpath = '2013/03/05/20140305_013000_59.fit.gz'
    d = datetime.datetime(2014, 3, 5, 1, 30)
    assert not s.matches(testpath, d)


def testDirectoryObsPattern():
    s = Scraper('%y%m%d/{observatory}_%Y%m%d.fits', observatory='SDO')
    testpath = '140305/SDO_20140305.fits'
    d = datetime.datetime(2014, 3, 5)
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
    startdate = datetime.datetime(2010, 10, 10, 5, 0)
    enddate = datetime.datetime(2010, 10, 10, 7, 0)
    timerange = TimeRange(startdate, enddate)
    assert len(s.range(timerange)) == 1


def testDirectoryRange_Month():
    s = Scraper('%Y%m/%d/%j_%H.txt')
    startdate = datetime.datetime(2008, 2, 20, 10)
    enddate = datetime.datetime(2008, 3, 2, 5)
    timerange = TimeRange(startdate, enddate)
    assert len(s.range(timerange)) == 12
    startdate = datetime.datetime(2009, 2, 20, 10)
    enddate = datetime.datetime(2009, 3, 2, 5)
    timerange = TimeRange(startdate, enddate)
    assert len(s.range(timerange)) == 11


def testNoDirectory():
    s = Scraper('files/%Y%m%d_%H%M.dat')
    startdate = datetime.datetime(2010, 1, 10, 20, 30)
    enddate = datetime.datetime(2010, 1, 20, 20, 30)
    timerange = TimeRange(startdate, enddate)
    assert len(s.range(timerange)) == 1


def testExtractDates_usingPattern():
    # Standard pattern
    s = Scraper('data/%Y/%m/%d/fits/swap/swap_00174_fd_%Y%m%d_%H%M%S.fts.gz')
    testURL = 'data/2014/05/14/fits/swap/swap_00174_fd_20140514_200135.fts.gz'
    timeURL = datetime.datetime(2014, 5, 14, 20, 1, 35)
    assert s._extractDateURL(testURL) == timeURL
    # Not-full repeated pattern
    s = Scraper('data/%Y/fits/swap/swap_00174_fd_%Y%m%d_%H%M%S.fts.gz')
    testURL = 'data/2014/fits/swap/swap_00174_fd_20140514_200135.fts.gz'
    timeURL = datetime.datetime(2014, 5, 14, 20, 1, 35)
    assert s._extractDateURL(testURL) == timeURL


def testExtractDates_notSeparators():
    s = Scraper('data/%Y/%m/swap%m%d_%H%M%S')
    testURL = 'data/2014/05/swap0514_200135'
    timeURL = datetime.datetime(2014, 5, 14, 20, 1, 35)
    assert s._extractDateURL(testURL) == timeURL


def testExtractDates_notSeparators_andSimilar():
    s = Scraper('data/%Y/Jun%b%d_%H%M%S')
    testURL = 'data/2014/JunJun14_200135'
    timeURL = datetime.datetime(2014, 6, 14, 20, 1, 35)
    assert s._extractDateURL(testURL) == timeURL
    testURL = 'data/2014/JunMay14_200135'
    timeURL = datetime.datetime(2014, 5, 14, 20, 1, 35)
    assert s._extractDateURL(testURL) == timeURL
    # and testing with the month afterwards
    s = Scraper('data/%Y/%dJun%b_%H%M%S')
    testURL = 'data/2014/14JunJun_200135'
    timeURL = datetime.datetime(2014, 6, 14, 20, 1, 35)
    assert s._extractDateURL(testURL) == timeURL


def testURL_pattern():
    s = Scraper('fd_%Y%m%d_%H%M%S.fts')
    assert s._URL_followsPattern('fd_20130410_231211.fts')
    assert not s._URL_followsPattern('fd_20130410_231211.fts.gz')
    assert not s._URL_followsPattern('fd_20130410_ar_231211.fts.gz')


@pytest.mark.xfail
def testURL_patternMilliseconds():
    s = Scraper('fd_%Y%m%d_%H%M%S_%e.fts')
    # NOTE: Seems that if below fails randomly - not understood why
    #       with `== True` fails a bit less...
    assert s._URL_followsPattern('fd_20130410_231211_119.fts')
    assert not s._URL_followsPattern('fd_20130410_231211.fts.gz')
    assert not s._URL_followsPattern('fd_20130410_ar_231211.fts.gz')

# Local files don't work
# def testFilesRange_sameDirectory_local():
#     s = Scraper('/'.join(['file:/',sunpy.data.test.rootdir,
#                           'EIT','efz%Y%m%d.%H%M%S_s.fits']))
#     print(s.pattern)
#     print(s.now)
#     startdate = datetime.datetime(2004, 3, 1, 4, 0)
#     enddate = datetime.datetime(2004, 3, 1, 6, 30)
#     assert len(s.filelist(TimeRange(startdate, enddate))) == 3
#     startdate = datetime.datetime(2010, 1, 10, 20, 30)
#     enddate = datetime.datetime(2010, 1, 20, 20, 30)
#     assert len(s.filelist(TimeRange(startdate, enddate))) == 0


@pytest.mark.xfail
@pytest.mark.remote_data
def testFilesRange_sameDirectory_remote():
    pattern = ('http://solarmonitor.org/data/%Y/%m/%d/'
               'fits/{instrument}/'
               '{instrument}_00174_fd_%Y%m%d_%H%M%S.fts.gz')
    s = Scraper(pattern, instrument='swap')
    startdate = datetime.datetime(2014, 5, 14, 0, 0)
    enddate = datetime.datetime(2014, 5, 14, 6, 30)
    timerange = TimeRange(startdate, enddate)
    assert len(s.filelist(timerange)) == 2
    startdate = datetime.datetime(2014, 5, 14, 21, 0)
    enddate = datetime.datetime(2014, 5, 14, 23, 30)
    timerange = TimeRange(startdate, enddate)
    assert len(s.filelist(timerange)) == 0


@pytest.mark.xfail
@pytest.mark.remote_data
def testFilesRange_sameDirectory_months_remote():
    pattern = ('http://www.srl.caltech.edu/{spacecraft}/DATA/{instrument}/'
               'Ahead/1minute/AeH%y%b.1m')
    s = Scraper(pattern, spacecraft='STEREO', instrument='HET')
    startdate = datetime.datetime(2007, 8, 1)
    enddate = datetime.datetime(2007, 9, 10)
    timerange = TimeRange(startdate, enddate)
    assert len(s.filelist(timerange)) == 2


@pytest.mark.remote_data
def test_ftp():
    pattern = 'ftp://solar-pub.nao.ac.jp/pub/nsro/norh/data/tcx/%Y/%m/tca%y%m%d'
    s = Scraper(pattern)
    timerange = TimeRange('2016/5/18 15:28:00', '2016/5/20 16:30:50')
    assert len(s.filelist(timerange)) == 2
