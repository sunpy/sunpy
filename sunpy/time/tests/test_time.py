from __future__ import absolute_import, division, print_function

from datetime import datetime

from sunpy import time
from sunpy.time import parse_time

import numpy as np
import pandas
from sunpy.extern.six.moves import range

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
    assert parse_time(765548612.0,'utime') == datetime(2003, 4, 5, 12, 23, 32)
    assert parse_time(1009685652.0,'utime') == datetime(2010, 12, 30, 4, 14, 12)

def test_parse_time_pandas_timestamp():
    ts = pandas.Timestamp(LANDING)

    dt = parse_time(ts)

    assert isinstance(dt, datetime)
    assert dt == LANDING

def test_parse_time_pandas_index():
    inputs = [datetime(2012, 1, i) for i in range(1,13)]
    ind = pandas.tseries.index.DatetimeIndex(inputs)

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
        ('20070504_210812', datetime(2007, 5, 4, 21, 8, 12))
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

