import logging
import datetime
from urllib.error import URLError, HTTPError
from unittest.mock import Mock, patch

import pytest

from sunpy.data.test import rootdir
from sunpy.extern import parse
from sunpy.net.scraper import Scraper
from sunpy.net.scraper_utils import get_timerange_from_exdict
from sunpy.time import TimeRange, parse_time


def test_directory_date_pattern():
    s = Scraper(format='{{year:4d}}/{{month:2d}}/{{day:2d}}/{{year:4d}}{{month:2d}}{{day:2d}}_{{hour:2d}}{{minute:2d}}{{second:2d}}_59.fit.gz')
    testpath = '2014/03/05/20140305_013000_59.fit.gz'
    d = parse_time((2014, 3, 5, 1, 30))
    assert s.matches(testpath, d)


def test_directory_date_patternFalse():
    s = Scraper(format='{{year:4d}}/{{month:2d}}/{{day:2d}}/{{year:4d}}{{month:2d}}{{day:2d}}_{{hour:2d}}{{minute:2d}}{{second:2d}}_59.fit.gz')
    testpath = '2013/03/05/20140305_013000_59.fit.gz'
    d = parse_time((2014, 3, 5, 1, 30))
    assert not s.matches(testpath, d)


def test_directory_obs_pattern():
    s = Scraper(format='{{year:2d}}{{month:2d}}{{day:2d}}/{observatory}_{{year:4d}}{{month:2d}}{{day:2d}}.fits', observatory='SDO')
    testpath = '140305/SDO_20140305.fits'
    d = parse_time((2014, 3, 5))
    assert s.matches(testpath, d)

def test_directory_range():
    s = Scraper(format='{{year:4d}}/{{month:2d}}/{{day:2d}}/{{year:4d}}{{month:2d}}{{day:2d}}_{{hour:2d}}.fit.gz')
    directory_list = ['2009/12/30/', '2009/12/31/', '2010/01/01/',
                      '2010/01/02/', '2010/01/03/']
    timerange = TimeRange('2009-12-30', '2010-01-03')
    assert s.range(timerange) == directory_list


def test_directory_regex():
    # Test for Windows where '\' is a path separator and not part of the regex
    s = Scraper(format='scheme://a.url.with/a/few/forward/slashes/andbacklash\\inthename.ext')
    timerange = TimeRange('2019-02-01', '2019-02-03')
    directory = s.range(timerange)
    assert directory == ['scheme://a.url.with/a/few/forward/slashes/']


def test_directory_range_false():
    s = Scraper(format='{{year:4d}}{{month:2d}}{{day:2d}}/{{year:4d}}{{month:2d}}{{day:2d}}_{{hour:2d}}.fit.gz')
    directory_list = ['20091230/', '20091231/', '20100101/',
                      '20090102/', '20090103/']
    timerange = TimeRange('2009/12/30', '2010/01/03')
    assert s.range(timerange) != directory_list


def testNoDateDirectory():
    s = Scraper(format='mySpacecraft/myInstrument/xMinutes/aaa{{year:4d}}{{month_name_abbr:l}}.ext')
    directory_list = ['mySpacecraft/myInstrument/xMinutes/']
    timerange = TimeRange('2009/11/20', '2010/01/03')
    assert s.range(timerange) == directory_list


def test_directory_range_hours():
    s = Scraper(format='{{year:4d}}{{month:2d}}{{day:2d}}_{{hour:2d}}/{{hour:2d}}{{minute:2d}}.csv')
    timerange = TimeRange('2009-12-31T23:40:00', '2010-01-01T01:15:00')
    assert len(s.range(timerange)) == 3  # 3 directories (1 per hour)


def test_directory_range_single():
    s = Scraper(format='{{year:4d}}{{month:2d}}{{day:2d}}/{{hour:2d}}_{{minute:2d}}.csv')
    startdate = parse_time((2010, 10, 10, 5, 0))
    enddate = parse_time((2010, 10, 10, 7, 0))
    timerange = TimeRange(startdate, enddate)
    assert len(s.range(timerange)) == 1


def test_directory_range_month():
    s = Scraper(format='{{year:4d}}{{month:2d}}/{{day:2d}}/{{day_of_year:3d}}_{{hour:2d}}.txt')
    startdate = parse_time((2008, 2, 20, 10))
    enddate = parse_time((2008, 3, 2, 5))
    timerange = TimeRange(startdate, enddate)
    assert len(s.range(timerange)) == 12
    startdate = parse_time((2009, 2, 20, 10))
    enddate = parse_time((2009, 3, 2, 5))
    timerange = TimeRange(startdate, enddate)
    assert len(s.range(timerange)) == 11


@pytest.mark.parametrize(('pattern', 'filename', 'metadict'), [
    ('_{{year:4d}}{{month:2d}}{{day:2d}}__{{millisecond:3d}}c{{:5d}}_{{:2d}}{{}}.fts', '_20201535__012c12345_33 .fts',
     {'year': 2020, 'month': 15, 'day': 35, 'millisecond': 12}),
    ('fd_{{year:4d}}{{month:2d}}{{day:2d}}_{{hour:2d}}{{minute:2d}}{{second:2d}}_{{millisecond:3d}}.fts', 'fd_20130410_231211_119.fts',
     {'year': 2013, 'month': 4, 'day': 10, 'hour': 23, 'minute': 12, 'second': 11, 'millisecond': 119})
])
def test_url_pattern(pattern, filename, metadict):
    assert parse(pattern.format(None), filename).named == metadict


def test_url_pattern_milliseconds_zero_padded():
    # Asserts solution to ticket #1954.
    # Milliseconds must be zero-padded in order to match URL lengths.
    now_mock = Mock(return_value=datetime.datetime(2019, 4, 19, 0, 0, 0, 4009))
    with patch('sunpy.net.scraper.datetime', now=now_mock):
        s = Scraper(format='fd_{{year:4d}}{{month:2d}}{{day:2d}}_{{hour:2d}}{{minute:2d}}{{second:2d}}_{{millisecond:3d}}.fts')
    now_mock.assert_called_once()
    assert s.now == 'fd_20190419_000000_004.fts'


def test_files_range_same_directory_local():
    s = Scraper(format='/'.join(['file:/', str(rootdir),
                          'EIT_header', 'efz{{year:4d}}{{month:2d}}{{day:2d}}.{{hour:2d}}{{minute:2d}}{{second:2d}}_s.header']))
    startdate = parse_time((2004, 3, 1, 4, 0))
    enddate = parse_time((2004, 3, 1, 6, 30))
    assert len(s.filelist(TimeRange(startdate, enddate))) == 3
    startdate = parse_time((2010, 1, 10, 20, 30))
    enddate = parse_time((2010, 1, 20, 20, 30))
    assert len(s.filelist(TimeRange(startdate, enddate))) == 0


@pytest.mark.remote_data
def test_files_range_same_directory_remote():
    pattern = ('https://proba2.sidc.be/{instrument}/data/bsd/{{year:4d}}/{{month:2d}}/{{day:2d}}/'
               '{instrument}_lv1_{{year:4d}}{{month:2d}}{{day:2d}}_{{hour:2d}}{{minute:2d}}{{second:2d}}.fits')
    s = Scraper(format=pattern, instrument='swap')
    startdate = parse_time((2014, 5, 14, 0, 0))
    enddate = parse_time((2014, 5, 14, 0, 5))
    timerange = TimeRange(startdate, enddate)
    assert len(s.filelist(timerange)) == 2
    startdate = parse_time((2014, 5, 14, 0, 6))
    enddate = parse_time((2014, 5, 14, 0, 7))
    timerange = TimeRange(startdate, enddate)
    assert len(s.filelist(timerange)) == 0


@pytest.mark.remote_data
def test_files_range_same_directory_months_remote():
    pattern = ('https://izw1.caltech.edu/{spacecraft}/DATA/{instrument}/'
               'Ahead/1minute/AeH{{year:2d}}{{month_name_abbr:w}}.1m')
    s = Scraper(format=pattern, spacecraft='STEREO', instrument='HET')
    startdate = parse_time((2007, 8, 1))
    enddate = parse_time((2007, 9, 10))
    timerange = TimeRange(startdate, enddate)
    files = s.filelist(timerange)
    assert files == ['https://izw1.caltech.edu/STEREO/DATA/HET/Ahead/1minute/AeH07Aug.1m',
                     'https://izw1.caltech.edu/STEREO/DATA/HET/Ahead/1minute/AeH07Sep.1m']


@pytest.mark.remote_data
def test_ftp():
    pattern = 'ftp://ftp.ngdc.noaa.gov/STP/space-weather/swpc-products/daily_reports/solar_region_summaries/{{year:4d}}/{{month:2d}}/{{year:4d}}{{month:2d}}{{day:2d}}SRS.txt'
    s = Scraper(format=pattern)
    timerange = TimeRange('2024/5/18', '2024/5/20')
    urls = s.filelist(timerange)
    assert urls[0] == ('ftp://ftp.ngdc.noaa.gov/STP/space-weather/swpc-products/daily_reports/solar_region_summaries/2024/05/20240518SRS.txt')
    assert len(urls) == 3


@pytest.mark.remote_data
def test_filelist_url_missing_directory():
    # Asserts solution to ticket #2684.
    # Attempting to access data for the year 1960 results in a 404, so no files are returned.
    pattern = 'https://lasp.colorado.edu/eve/data_access/evewebdata/products/level2/{{year:4d}}/{{day_of_year:3d}}/'
    s = Scraper(format=pattern)
    timerange = TimeRange('1960/01/01 00:00:00', '1960/01/02 00:00:00')
    assert len(s.filelist(timerange)) == 0


@pytest.mark.remote_data
def test_filelist_relative_hrefs():
    # the url opened by the scraper from below pattern contains some links which don't have hrefs
    pattern = 'https://www.bbso.njit.edu/pub/archive/{{year:4d}}/{{month:2d}}/{{day:2d}}/bbso_halph_fr_{{year:4d}}{{month:2d}}{{day:2d}}_{{hour:2d}}{{minute:2d}}{{second:2d}}.fts'
    s = Scraper(format=pattern)
    timerange = TimeRange('2016/5/18 15:28:00', '2016/5/18 16:30:00')
    assert s.domain == 'https://www.bbso.njit.edu/'
    # hrefs are relative to domain here, not to the directory they are present in
    # this checks that `scraper.filelist` returns fileurls relative to the domain
    fileurls = s.filelist(timerange)
    assert fileurls[1] == s.domain + 'pub/archive/2016/05/18/bbso_halph_fr_20160518_160033.fts'


@pytest.mark.remote_data
def test_extract_files_meta():
    pattern0 = 'https://solar.nro.nao.ac.jp/norh/data/tcx/{{year:4d}}/{{month:2d}}/{{wave}}{{year:2d}}{{month:2d}}{{day:2d}}'
    s0 = Scraper(format=pattern0)
    timerange0 = TimeRange('2020/1/1 4:00', '2020/1/2')
    matchdict = {'wave': ['tca', 'tcz']}
    metalist0 = s0._extract_files_meta(timerange0, matcher=matchdict)
    assert metalist0[0]['wave'] == 'tca'
    assert metalist0[3]['wave'] == 'tcz'
    assert metalist0[1]['day'] == 2

    pattern1 = ('https://gong2.nso.edu/oQR/zqs/{{year:4d}}{{month:2d}}/mrzqs{{year:2d}}{{month:2d}}{{day:2d}}'
                '/mrzqs{{year:2d}}{{month:2d}}{{day:2d}}t{{hour:2d}}{{minute:2d}}c{{CAR_ROT:4d}}_{{:3d}}.fits.gz')
    s1 = Scraper(format=pattern1)
    timerange1 = TimeRange('2020-01-05', '2020-01-05T16:00:00')
    metalist1 = s1._extract_files_meta(timerange1)
    urls = s1.filelist(timerange1)
    assert metalist1[3]['CAR_ROT'] == 2226
    assert metalist1[-1]['url'] == urls[-1]


def test_no_directory():
    s = Scraper(format='files/{{year:4d}}{{month:2d}}{{day:2d}}_{{hour:2d}}{{month:2d}}.dat')
    startdate = parse_time((2010, 1, 10, 20, 30))
    enddate = parse_time((2010, 1, 20, 20, 30))
    timerange = TimeRange(startdate, enddate)
    assert len(s.range(timerange)) == 1


@pytest.mark.parametrize(('exdict', 'start', 'end'), [
    ({"year": 2000}, '2000-01-01 00:00:00', '2000-12-31 23:59:59.999000'),
    ({"year": 2016, "month": 2}, '2016-02-01 00:00:00', '2016-02-29 23:59:59.999000'),
    ({'year': 2019, 'month': 2, 'day': 28}, '2019-02-28 00:00:00', '2019-02-28 23:59:59.999000'),
    ({'year': 2020, 'month': 7, 'day': 31, 'hour': 23, 'minute': 59, 'second': 59},
     '2020-07-31 23:59:59', '2020-07-31 23:59:59.999000')])
def test_get_timerange_with_extractor(exdict, start, end):
    tr = TimeRange(start, end)
    file_timerange = get_timerange_from_exdict(exdict)
    assert file_timerange == tr


@pytest.mark.remote_data
def test_parse_pattern_data():
    prefix = 'https://gong2.nso.edu/oQR/zqs/'
    pattern = prefix + '{{year:4d}}{{month:2d}}/mrzqs{{year:2d}}{{month:2d}}{{day:2d}}/mrzqs{{year:2d}}{{month:2d}}{{day:2d}}t{{hour:2d}}{{minute:2d}}c{{:4d}}_{{:3d}}.fits.gz'
    s = Scraper(format=pattern)
    timerange = TimeRange('2020-01-05', '2020-01-06T16:00:00')
    assert len(s.filelist(timerange)) == 37


@pytest.mark.remote_data
def test_yearly_overlap():
    # Check that a time range that falls within the interval that a file represents
    # returns a single result.
    pattern = "https://www.ngdc.noaa.gov/stp/space-weather/solar-data/solar-features/solar-flares/x-rays/goes/xrs/goes-xrs-report_{{year:4d}}.txt"
    scraper = Scraper(format=pattern)

    # Should return a single file for 2013
    trange = TimeRange("2013-01-02", "2013-01-03")
    assert len(scraper.filelist(trange)) == 1


@pytest.mark.parametrize(
    ("error_code", "expected_number_calls", "error_message"),
    [
        (400, 1, "Bad Request"),
        (403, 1, "Forbidden"),
        (429, 6, "Too Many Requests"),
        (504, 6, "Gateway Timeout"),
    ],
)
def test_http_errors_with_enqueue_limit(error_code, expected_number_calls, error_message):
    with patch("sunpy.net.scraper.urlopen") as mocked_urlopen:
        mocked_urlopen.side_effect = HTTPError(
            "https://example.com", error_code, error_message, {}, None
        )
        time_range = TimeRange("2012/3/4", "2012/3/4 02:00")
        pattern = "https://proba2.sidc.be/lyra/data/bsd/{{year:4d}}/{{month:2d}}/{{day:2d}}/{{}}_lev{{Level:1d}}_std.fits"
        scraper = Scraper(format=pattern)
        with pytest.raises(HTTPError, match=error_message) as excinfo:
            scraper._httpfilelist(time_range)
        assert excinfo.value.code == error_code
        assert mocked_urlopen.call_count == expected_number_calls


def test_connection_error():
    with patch('sunpy.net.scraper.urlopen') as mocked_urlopen:
        mocked_urlopen.side_effect = URLError('connection error')
        time = TimeRange('2012/3/4', '2012/3/4 02:00')
        pattern = "https://proba2.sidc.be/lyra/data/bsd/{{year:4d}}/{{month:2d}}/{{day:2d}}/{{}}_lev{{Level:1d}}_std.fits"
        scraper = Scraper(format=pattern)
        with pytest.raises(URLError, match='connection error'):
            scraper._httpfilelist(time)


def test_http_404_error_debug_message(caplog):
    with caplog.at_level(logging.DEBUG, logger='sunpy'):
        def patch_range(self, range):
            return ['https://test.com/']
        with patch('sunpy.net.scraper.urlopen') as mocked_urlopen:
            with patch.object(Scraper, 'range', patch_range):
                mocked_urlopen.side_effect = HTTPError('https://example.com', 404, '', {}, None)
                time = TimeRange('2012/3/4', '2012/3/4 02:00')
                pattern = "https://proba2.sidc.be/lyra/data/bsd/{{year:4d}}/{{month:2d}}/{{day:2d}}/{{}}_lev{{Level:1d}}_std.fits"
                scraper = Scraper(format=pattern)
                scraper._httpfilelist(time)
                assert "Directory https://test.com/ not found." in caplog.text


def test_check_timerange():
    s = Scraper(format='{{year:4d}}.fits')
    # Valid time range for 2014.fits is the whole of 2014
    # Test different cases to make sure check_timerange is working as expected

    # Interval exactly on lower boundary
    assert s._check_timerange('2014.fits', TimeRange("2013-06-01", "2014-01-01"))
    # Overlaps lower boundary
    assert s._check_timerange('2014.fits', TimeRange("2013-06-01", "2014-01-02"))
    # Overlaps upper and lower boundary
    assert s._check_timerange('2014.fits', TimeRange("2013-06-01", "2015-01-02"))
    # Entirely within both boundaries
    assert s._check_timerange('2014.fits', TimeRange("2014-06-01", "2014-07-02"))
    # Overlaps upper boundary
    assert s._check_timerange('2014.fits', TimeRange("2014-06-01", "2015-01-02"))
    # Interval exactly on upper boundary
    assert not s._check_timerange('2014.fits', TimeRange("2015-01-01", "2015-01-02"))

    # Interval below both boundaries
    assert not s._check_timerange('2014.fits', TimeRange("2002-01-01", "2013-01-02"))
    # Interval above both boundaries
    assert not s._check_timerange('2014.fits', TimeRange("2022-01-01", "2025-01-02"))

    # Only 2-digit Year
    assert s._check_timerange('14.fits', TimeRange("2013-06-01", "2014-01-01"))


def test_local_expected_directory_doesnt_exist(tmp_path):
    path = (tmp_path / '2025' / '01' / '01')
    path.mkdir(parents=True)
    with (path / 'test.txt').open('w') as file:
        file.write('')
    s = Scraper(format='file://'+str(tmp_path)+'/{{year:4d}}/{{month:2d}}/{{day:2d}}/{{file}}')
    files = s.filelist(TimeRange("2025-01-01", "2025-01-02"))
    assert len(files) == 1
