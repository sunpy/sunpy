import pytest
from dateutil.relativedelta import relativedelta

from sunpy.net.scraper_utils import check_timerange, get_timerange_from_exdict, smaller_pattern
from sunpy.time import TimeRange

TIMEPATTERN_EXAMPLES = [
    ('%b%y', relativedelta(months=1)),
    ('%m%y', relativedelta(months=1)),
    ('%H%d', relativedelta(hours=1)),
    ('%y%b', relativedelta(months=1)),
]

@pytest.mark.parametrize(('pattern', 'mintime'), TIMEPATTERN_EXAMPLES)
def test_smaller_pattern(pattern, mintime):
    assert mintime == smaller_pattern(pattern)


def test_check_timerange():
    pattern = '{{year:4d}}.fits'
    # Valid time range for 2014.fits is the whole of 2014
    # Test different cases to make sure check_timerange is working as expected

    # Interval exactly on lower boundary
    assert check_timerange(pattern, '2014.fits', TimeRange("2013-06-01", "2014-01-01"))
    # Overlaps lower boundary
    assert check_timerange(pattern, '2014.fits', TimeRange("2013-06-01", "2014-01-02"))
    # Overlaps upper and lower boundary
    assert check_timerange(pattern, '2014.fits', TimeRange("2013-06-01", "2015-01-02"))
    # Entirely within both boundaries
    assert check_timerange(pattern, '2014.fits', TimeRange("2014-06-01", "2014-07-02"))
    # Overlaps upper boundary
    assert check_timerange(pattern, '2014.fits', TimeRange("2014-06-01", "2015-01-02"))
    # Interval exactly on upper boundary
    assert not check_timerange(pattern, '2014.fits', TimeRange("2015-01-01", "2015-01-02"))

    # Interval below both boundaries
    assert not check_timerange(pattern, '2014.fits', TimeRange("2002-01-01", "2013-01-02"))
    # Interval above both boundaries
    assert not check_timerange(pattern, '2014.fits', TimeRange("2022-01-01", "2025-01-02"))

    pattern_month_name = '{{year:4d}}-{{month_name:l}}.fits'
    assert check_timerange(pattern_month_name, '2014-March.fits', TimeRange("2014-03-01", "2004-04-01"))

    pattern_month_name_abbr = '{{year:4d}}-{{month_name_abbr:l}}-{{day:2d}}.fits'
    assert check_timerange(pattern_month_name_abbr, '2004-Mar-06.fits', TimeRange("2004-03-06", "2004-03-07"))

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
