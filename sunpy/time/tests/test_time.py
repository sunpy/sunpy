from __future__ import absolute_import, division, print_function

from datetime import datetime

from sunpy import time
from sunpy.time import parse_time, is_time_in_given_format, get_day, find_time

import astropy.time
import numpy as np
import pandas
from sunpy.extern.six.moves import range

import pytest

LANDING = datetime(1966, 2, 3)


def test_parse_time_24():
    assert parse_time("2010-10-10T24:00:00") == datetime(2010, 10, 11)


def test_parse_time_24_2():
    assert parse_time("2010-10-10T24:00:00.000000") == datetime(2010, 10, 11)


def test_parse_time_trailing_zeros():
    # see issue #289 at https://github.com/sunpy/sunpy/issues/289
    assert parse_time('2010-10-10T00:00:00.00000000') == datetime(2010, 10, 10)


def test_parse_time_tuple():
    assert parse_time((1966, 2, 3)) == LANDING


def test_parse_time_int():
    assert parse_time(765548612.0, 'utime') == datetime(2003, 4, 5,
                                                        12, 23, 32)
    assert parse_time(1009685652.0, 'utime') == datetime(2010, 12, 30,
                                                         4, 14, 12)


def test_parse_time_pandas_timestamp():
    ts = pandas.Timestamp(LANDING)

    dt = parse_time(ts)

    assert isinstance(dt, datetime)
    assert dt == LANDING


def test_parse_time_pandas_series():
    inputs = [datetime(2012, 1, i) for i in range(1, 13)]
    ind = pandas.Series(inputs)

    dts = parse_time(ind)

    assert isinstance(dts, np.ndarray)
    assert all([isinstance(dt, datetime) for dt in dts])
    assert list(dts) == inputs


def test_parse_time_pandas_index():
    inputs = [datetime(2012, 1, i) for i in range(1, 13)]
    ind = pandas.DatetimeIndex(inputs)

    dts = parse_time(ind)

    assert isinstance(dts, np.ndarray)
    assert all([isinstance(dt, datetime) for dt in dts])
    assert list(dts) == inputs


def test_parse_time_numpy_date():
    inputs = np.arange('2005-02', '2005-03', dtype='datetime64[D]')

    dts = parse_time(inputs)

    assert isinstance(dts, np.ndarray)
    assert all([isinstance(dt, datetime) for dt in dts])


def test_parse_time_numpy_datetime():
    inputs = np.arange('2005-02-01T00', '2005-02-01T10', dtype='datetime64')

    dts = parse_time(inputs)

    assert isinstance(dts, np.ndarray)
    assert all([isinstance(dt, datetime) for dt in dts])


def test_parse_time_individual_numpy_datetime():
    dt64 = np.datetime64('2005-02-01T00')
    dt = parse_time(dt64)

    assert isinstance(dt, datetime)


def test_parse_time_numpy_datetime_timezone():
    dt64 = np.datetime64('2014-02-07T16:47:51-0500')
    dt = parse_time(dt64)

    assert dt == datetime(2014, 2, 7, 21, 47, 51)


def test_parse_time_numpy_datetime_ns():
    dt64 = np.datetime64('2014-02-07T16:47:51.008288000')
    dt = parse_time(dt64)

    assert dt == datetime(2014, 2, 7, 16, 47, 51, 8288)

    dt64 = np.datetime64('2014-02-07T16:47:51.008288123')
    dt = parse_time(dt64)

    assert dt == datetime(2014, 2, 7, 16, 47, 51, 8288)


def test_parse_time_numpy_datetime_round():
    dt64 = np.datetime64('2014-02-07T16:47:51.008288999')
    dt = parse_time(dt64)

    assert dt == datetime(2014, 2, 7, 16, 47, 51, 8288)


def test_parse_time_astropy():
    astropy_time = parse_time(astropy.time.Time(['2016-01-02T23:00:01']))

    assert astropy_time == datetime(year=2016, month=1, day=2, hour=23, minute=0, second=1)


def test_parse_time_now():
    """
    Ensure 'parse_time' can be called with 'now' argument to get utc
    """
    # TODO: once mocking support is merged in, we can perform a like for like comparison,
    #       the following at least ensures that 'now' is a legal argument.
    assert isinstance(parse_time('now'), datetime) is True


def test_ISO():
    assert parse_time('1966-02-03') == LANDING
    assert (
        parse_time('1966-02-03T20:17:40') == datetime(1966, 2, 3, 20, 17, 40)
    )
    assert (
        parse_time('19660203T201740') == datetime(1966, 2, 3, 20, 17, 40)
    )

    lst = [
        ('2007-05-04T21:08:12.999999',
         datetime(2007, 5, 4, 21, 8, 12, 999999)),
        ('20070504T210812.999999',
         datetime(2007, 5, 4, 21, 8, 12, 999999)),
        ('2007/05/04 21:08:12.999999',
         datetime(2007, 5, 4, 21, 8, 12, 999999)),
        ('2007-05-04 21:08:12.999999',
         datetime(2007, 5, 4, 21, 8, 12, 999999)),
        ('2007/05/04 21:08:12', datetime(2007, 5, 4, 21, 8, 12)),
        ('2007-05-04 21:08:12', datetime(2007, 5, 4, 21, 8, 12)),
        ('2007-05-04 21:08', datetime(2007, 5, 4, 21, 8)),
        ('2007-05-04T21:08:12', datetime(2007, 5, 4, 21, 8, 12)),
        ('20070504T210812', datetime(2007, 5, 4, 21, 8, 12)),
        ('2007-May-04 21:08:12', datetime(2007, 5, 4, 21, 8, 12)),
        ('2007-May-04 21:08', datetime(2007, 5, 4, 21, 8)),
        ('2007-May-04', datetime(2007, 5, 4)),
        ('2007-05-04', datetime(2007, 5, 4)),
        ('2007/05/04', datetime(2007, 5, 4)),
        ('04-May-2007', datetime(2007, 5, 4)),
        ('04-May-2007 21:08:12.999999', datetime(2007, 5, 4, 21, 8, 12, 999999)),
        ('20070504_210812', datetime(2007, 5, 4, 21, 8, 12)),
        ('2007.05.04_21:08:12_TAI', datetime(2007, 5, 4, 21, 8, 12)),
    ]

    for k, v in lst:
        assert parse_time(k) == v


def test_break_time():
    t = datetime(2007, 5, 4, 21, 8, 12)
    assert time.break_time(t) == '20070504_210812'


def test_day_of_year():
    # Note that 2012 is a leap year, 2011 is a standard year
    # test that it starts at 1
    assert time.day_of_year('2011/01/01') == 1.0
    # test fractional day
    assert time.day_of_year('2011/01/01 06:00') == 1.25
    assert time.day_of_year('2011/01/01 12:00') == 1.50
    assert time.day_of_year('2011/01/01 18:00') == 1.75
    # test correct number of days in a (standard) year
    assert time.day_of_year('2011/12/31') == 365
    # test correct number of days in a (leap) year
    assert time.day_of_year('2012/12/31') == 366
    # test a few extra dates in standard year
    assert time.day_of_year('2011/08/01') == 213
    assert time.day_of_year('2011/04/10') == 100
    assert time.day_of_year('2011/01/31') == 31
    assert time.day_of_year('2011/09/30') == 273
    # test a few extra dates in a leap year
    assert time.day_of_year('2012/08/01') == 214
    assert time.day_of_year('2012/04/10') == 101
    assert time.day_of_year('2012/01/31') == 31
    assert time.day_of_year('2012/09/30') == 274


def test_time_string_parse_format():
    assert parse_time('01/06/2012',
                      _time_string_parse_format='%d/%m/%Y') == datetime(2012, 6, 1, 0, 0)
    assert parse_time('06/01/2012',
                      _time_string_parse_format='%d/%m/%Y') == datetime(2012, 1, 6, 0, 0)
    assert parse_time('06/01/85',
                      _time_string_parse_format='%d/%m/%y') == datetime(1985, 1, 6, 0, 0)
    assert parse_time('6/1/85',
                      _time_string_parse_format='%d/%m/%y') == datetime(1985, 1, 6, 0, 0)
    with pytest.raises(ValueError):
        parse_time('01/06/2012')
    with pytest.raises(ValueError):
        parse_time('01/06/2012', time_string_parse_format='%d/%m/%m')

    with pytest.raises(ValueError):
        parse_time('2016', _time_string_parse_format='zz')


def test__iter_empty():

    class CountDown(object):

        def __init__(self, start_from=0):
            self.start = start_from

        def __iter__(self):
            return self

        def __next__(self):
            self.start -= 1

            if self.start < 0:
                raise StopIteration

            return self.start

        next = __next__   # Support Py2.x

    one_count = CountDown(1)
    assert time.time._iter_empty(one_count) is False
    assert time.time._iter_empty(one_count) is True



def test_is_time():
    assert time.is_time(datetime.utcnow()) is True
    assert time.is_time('2017-02-14 08:08:12.999', "%Y-%m-%d %H:%M:%S.%f") is True

    assert time.is_time(None) is False
    assert time.is_time('2016-14-14 19:08', "%Y-%b-%d %H:%M:%S") is False


def test_is_time_in_given_format():
    assert is_time_in_given_format('2017-02-14 08:08:12.999', "%Y-%m-%d %H:%M:%S.%f") is True
    assert is_time_in_given_format('2017-02-14 08:08:12.999', "%Y-%m-%dT%H:%M:%S.%f") is False


def test_get_day():
    end_of_day = datetime(year=2017, month=1, day=1, hour=23, minute=59, second=59,
                          microsecond=999)

    begining_of_day = get_day(end_of_day)
    assert begining_of_day.year == 2017
    assert begining_of_day.month == 1
    assert begining_of_day.day == 1
    assert begining_of_day.hour == 0
    assert begining_of_day.minute == 0
    assert begining_of_day.second == 0
    assert begining_of_day.microsecond == 0
