from datetime import UTC, date, datetime

import numpy as np
import pandas
import pytest

import astropy.time
from astropy.time import Time

import sunpy.time as time
from sunpy.time import is_time_equal, parse_time

LANDING = Time('1966-02-03', format='isot')


def test_parse_time_24():
    dt = parse_time("2010-10-10T24:00:00")
    assert dt == Time('2010-10-11')
    assert dt.format == 'isot'
    assert dt.scale == 'utc'


def test_parse_time_24_2():
    dt = parse_time("2010-10-10T24:00:00.000000")
    assert dt == Time('2010-10-11')
    assert dt.format == 'isot'
    assert dt.scale == 'utc'


def test_parse_time_microseconds_excess_trailing_zeros():
    dt = parse_time('2010-Oct-10 00:00:00.1234560')
    assert dt == Time('2010-10-10 00:00:00.123456')
    assert dt.format == 'isot'
    assert dt.scale == 'utc'

    # Excess digits beyond 6 digits should error if they are not zeros
    with pytest.raises(ValueError, match="Input values did not match any of the formats where the format keyword is optional:"):
        dt = parse_time('2010-Oct-10 00:00:00.1234567')

    # An ending run of zeros should still error if they are not a microsecond field
    with pytest.raises(ValueError, match="Input values did not match any of the formats where the format keyword is optional:"):
        dt = parse_time('10-Oct-2010.0000000')


def test_parse_time_tuple():
    dt = parse_time((1966, 2, 3))
    assert dt == LANDING
    assert dt.format == 'isot'
    assert dt.scale == 'utc'

    dt = parse_time((1966, 2, 3, 12, 2, 3))
    assert dt == Time('1966-2-3T12:2:3')
    assert dt.format == 'isot'
    assert dt.scale == 'utc'

    dt = parse_time((1966, 2, 3, 12, 2, 3, 8266))
    assert dt == Time('1966-2-3T12:2:3.008266')
    assert dt.format == 'isot'
    assert dt.scale == 'utc'


def test_parse_time_int():
    # Once https://github.com/astropy/astropy/issues/6970 is fixed,
    # remove .jd from equality check
    dt1 = parse_time(765548612.0, format='utime')
    assert dt1.jd == Time('2003-4-5T12:23:32').jd
    assert dt1.format == 'utime'

    dt2 = parse_time(1009685652.0, format='utime')
    assert dt2.jd == Time('2010-12-30T4:14:12').jd
    assert dt2.format == 'utime'


def test_parse_time_pandas_timestamp():
    ts = pandas.Timestamp(LANDING.datetime)

    dt = parse_time(ts)

    assert isinstance(dt, astropy.time.Time)
    assert dt == LANDING


def test_parse_time_nanoseconds():
    # Check that nanosecon precision is retained when parsing pandas timestamps
    ts = pandas.Timestamp('2020-07-31 00:00:26.166196864')
    dt = parse_time(ts)
    assert dt.jd1 == 2459062.0
    # If nanoseconds are not retained, this value is slightly too low
    assert dt.jd2 == -0.4996971504992593


def test_parse_time_pandas_series():
    inputs = [datetime(2012, 1, i) for i in range(1, 13)]
    ind = pandas.Series(inputs)
    as_inps = Time(inputs)

    dts = parse_time(ind)

    assert isinstance(dts, astropy.time.Time)
    assert np.all(dts == as_inps)


def test_parse_time_pandas_series_2():
    inputs = [[datetime(2012, 1, 1, 0, 0), datetime(2012, 1, 2, 0, 0)],
              [datetime(2012, 1, 3, 0, 0), datetime(2012, 1, 4, 0, 0)]]
    ind = pandas.Series(inputs)
    as_inps = Time(inputs)

    apts = parse_time(ind)

    assert isinstance(apts, astropy.time.Time)
    assert np.all(apts == as_inps)
    assert apts.shape == as_inps.shape


def test_parse_time_pandas_index():
    inputs = [datetime(2012, 1, i) for i in range(1, 13)]
    ind = pandas.DatetimeIndex(inputs)
    as_inps = Time(inputs)

    dts = parse_time(ind)

    assert isinstance(dts, astropy.time.Time)
    assert np.all(dts == as_inps)


def test_parse_time_numpy_date():
    inputs = np.arange('2005-02', '2005-03', dtype='datetime64[D]')

    dts = parse_time(inputs)

    assert isinstance(dts, astropy.time.Time)
    assert np.all(dts == Time([str(dt.astype('M8[ns]')) for dt in inputs]))


def test_parse_time_numpy_datetime():
    inputs = np.arange('2005-02-01T00', '2005-02-01T10', dtype='datetime64')

    dts = parse_time(inputs)

    assert isinstance(dts, astropy.time.Time)
    assert np.all(dts == Time([str(dt.astype('M8[ns]')) for dt in inputs]))


def test_parse_time_individual_numpy_datetime():
    dt64 = np.datetime64('2005-02-01T00')
    dt = parse_time(dt64)

    assert isinstance(dt, astropy.time.Time)
    assert dt == Time('2005-02-01', format='isot')


def test_parse_time_numpy_datetime_ns():
    dt64 = np.datetime64('2014-02-07T16:47:51.008288000')
    dt = parse_time(dt64)

    assert dt == Time('2014-02-07T16:47:51.008288000', format='isot')

    dt64 = np.datetime64('2014-02-07T16:47:51.008288123')
    dt = parse_time(dt64)

    assert dt == Time('2014-02-07T16:47:51.008288123', format='isot')

    dt64 = np.datetime64('2014-02-07T16:47:51.234565999')
    dt = parse_time(dt64)

    assert dt == Time('2014-02-07T16:47:51.234565999')


def test_parse_time_astropy():
    ip = astropy.time.Time(['2016-01-02T23:00:01'])
    astropy_time = parse_time(ip)

    assert astropy_time == ip
    assert astropy_time.format == 'isot'


def test_parse_time_datetime():
    dt = datetime(2014, 2, 7, 16, 47, 51, 8288)
    assert parse_time(dt) == Time('2014-02-07 16:47:51.008288')
    assert parse_time(dt).format == 'datetime'


def test_parse_time_date():
    dt = parse_time(date(1966, 2, 3))
    assert dt == Time('1966-2-3')
    assert dt.format == 'iso'


def test_parse_time_now():
    now = parse_time('now')
    assert isinstance(now, astropy.time.Time)
    assert now.format == 'datetime'
    assert now.scale == 'utc'


@pytest.mark.parametrize(('expected', 'time_string'), [
    ('2007-05-04T21:08:12.999999', '2007-05-04T21:08:12.999999'),
    ('2007-05-04T21:08:12.999999', '2007/05/04T21:08:12.999999'),
    ('2007-05-04T21:08:12.999999', '2007-05-04T21:08:12.999999Z'),
    ('2007-05-04T21:08:12.999999', '20070504T210812.999999'),
    ('2007-05-04T21:08:12.999999', '2007/05/04 21:08:12.999999'),
    ('2007-05-04T21:08:12.999999', '2007-05-04 21:08:12.999999'),
    ('2007-05-04T21:08:12.999999', '2007-May-04 21:08:12.999999'),
    ('2007-05-04T21:08:12.999999', '04-May-2007 21:08:12.999999'),
    ('2007-05-04T21:08:12.999999', '2007:124:21:08:12.999999'),
    ('2007-05-04T21:08:12', '2007-05-04T21:08:12'),
    ('2007-05-04T21:08:12', '2007/05/04T21:08:12'),
    ('2007-05-04T21:08:12', '2007-05-04T21:08:12Z'),
    ('2007-05-04T21:08:12', '20070504T210812'),
    ('2007-05-04T21:08:12', '2007/05/04 21:08:12'),
    ('2007-05-04T21:08:12', '2007-05-04 21:08:12'),
    ('2007-05-04T21:08:12', '2007-May-04 21:08:12'),
    ('2007-05-04T21:08:12', '04-May-2007 21:08:12'),
    ('2007-05-04T21:08:12', '20070504_210812'),
    ('2007-05-04T21:08:12', '2007:124:21:08:12'),
    ('2007-05-04T21:08:12', '20070504210812'),
    ('2007-05-04T21:08:12', '2007.05.04_21:08:12_UTC'),
    ('2007-05-04T21:08:12', '2007.05.04_21:08:12'),
    ('2007-05-04T21:08:00', '20070504T2108'),
    ('2007-05-04T21:08:00', '2007/05/04 21:08'),
    ('2007-05-04T21:08:00', '2007-05-04 21:08'),
    ('2007-05-04T21:08:00', '2007-May-04 21:08'),
    ('2007-05-04T21:08:00', '20070504_2108'),
    ('2007-05-04T21:08:00', '200705042108'),
    ('2007-05-04T21:08:00', '2007/05/04T21:08'),
    ('2007-05-04', '2007-May-04'),
    ('2007-05-04', '2007-05-04'),
    ('2007-05-04', '2007/05/04'),
    ('2007-05-04', '04-May-2007'),
])
def test_parse_time_utc(expected, time_string):
    dt = parse_time(time_string)
    assert is_time_equal(dt, Time(expected))
    assert dt.scale == 'utc'
    assert dt.format == 'isot'

    # Confirm that passing the string within a list gives the same answer
    assert is_time_equal(dt, parse_time([time_string])[0])


@pytest.mark.parametrize(('expected', 'time_string'), [
    ('2007-05-04T21:08:12.999999', '2007.05.04_21:08:12.999999_TAI'),
    ('2007-05-04T21:08:12', '2007.05.04_21:08:12_TAI'),
])
def test_parse_time_tai(expected, time_string):
    dt = parse_time(time_string)
    assert is_time_equal(dt, Time(expected, scale='tai'))
    assert dt.scale == 'tai'
    assert dt.format == 'isot'

    # Confirm that passing the string within a list gives the same answer
    assert is_time_equal(dt, parse_time([time_string])[0])


def test_parse_astropy_fallback():
    # If none of our string formats match, parse_time should simply pass through to Time
    assert parse_time('J2000') == Time('J2000')
    assert parse_time(['J2000']) == Time(['J2000'])


def test_parse_time_leap_second():
    dt1 = parse_time('1995-12-31 23:59:60')
    dt2 = Time('1995-12-31T23:59:60')

    assert dt1.jd == dt2.jd

    dt3 = parse_time('1995-Dec-31 23:59:60')

    assert dt2.jd == dt3.jd


# This warning shouldn't really be ignored; see https://github.com/astropy/astropy/issues/10564
# for the issue that will help decide how to handle this
@pytest.mark.filterwarnings('ignore:FITS time strings should no longer have embedded time scale')
@pytest.mark.parametrize(('ts', 'fmt'), [
    (1950.0, 'byear'),
    ('B1950.0', 'byear_str'),
    (63072064.184, 'cxcsec'),
    (datetime(2000, 1, 2, 12, 0, 0), 'datetime'),
    (2000.45, 'decimalyear'),
    ('2000-01-01T00:00:00.000(TAI)', 'fits'),
    (630720013.0, 'gps'),
    ('2000-01-01 00:00:00.000', 'iso'),
    ('2000-01-01T00:00:00.000', 'isot'),
    (2451544.5, 'jd'),
    (2000.0, 'jyear'),
    ('J2000.0', 'jyear_str'),
    (51544.0, 'mjd'),
    (730120.0003703703, 'plot_date'),
    (946684800.0, 'unix'),
    ('2000:001:00:00:00.000', 'yday')
])
def test_parse_time_astropy_formats(ts, fmt):
    dt = parse_time(ts, format=fmt)
    assert dt.format == fmt


def test_parse_time_int_float():
    # int and float values are not unique
    # The format has to be mentioned
    with pytest.raises(ValueError, match="Input values did not match any of the formats where the format keyword is optional:"):
        parse_time(100)
    with pytest.raises(ValueError, match="Input values did not match any of the formats where the format keyword is optional:"):
        parse_time(100.0)


@pytest.mark.parametrize("scale", [
    'tai',
    'tcb',
    'tcg',
    'tdb',
    'tt',
    'ut1',
    'utc'
])
def test_parse_time_scale(scale):
    dt = parse_time('2007-05-04T21:08:12', scale=scale)
    dt2 = Time('2007-05-04T21:08:12', scale=scale)
    assert is_time_equal(dt, dt2)
    assert dt.scale == scale

    dt = parse_time(np.datetime64('2007-05-04T21:08:12'), scale=scale)
    dt2 = Time('2007-05-04T21:08:12', scale=scale)
    assert dt == dt2
    assert dt.scale == scale

    dt = datetime(2014, 2, 7, 16, 47, 51)
    dt = parse_time(dt, scale=scale)
    dt2 = Time('2014-02-07T16:47:51', scale=scale)
    assert dt == dt2
    assert dt.scale == scale

    dt = date(2014, 2, 7)
    dt = parse_time(dt, scale=scale)
    dt2 = Time('2014-02-07', scale=scale)
    assert dt == dt2
    assert dt.scale == scale


def test_parse_time_list():
    tstrings = ['2010-09-03 00:00:00', '2005-09-03 06:00:00',
                '1995-12-31 23:59:60']

    assert np.all(parse_time(tstrings) == Time(tstrings))


def test_parse_time_list_2():
    tstrings = [['1998-01-01 00:00:01', '1998-01-01 00:00:02'],
                ['1998-01-01 00:00:03', '1995-12-31 23:59:60']]

    assert np.all(parse_time(tstrings) == Time(tstrings))


def test_parse_time_list_3():
    tstrings = ['2001-Jan-01', '2001-Jan-02', '2001-Jan-03']
    assert np.all(parse_time(tstrings) == Time.strptime(tstrings, '%Y-%b-%d'))


def test_parse_time_long_list():
    # Check whether lists of >500 strings cause a segmentation fault
    # https://github.com/astropy/astropy/issues/18254
    assert isinstance(parse_time(['2001-Jan-01']*501), Time)  # sunpy format
    assert isinstance(parse_time(['J2000.0']*501), Time)  # non-sunpy format


def test_is_time():
    time.is_time(datetime.now(UTC)) is True
    assert time.is_time('2017-02-14 08:08:12.999') is True
    assert time.is_time(Time.now()) is True

    assert time.is_time(None) is False
    assert time.is_time('2016-14-14 19:08') is False


def test_is_time_in_given_format():
    assert time.is_time_in_given_format('2017-02-14 08:08:12.999', "%Y-%m-%d %H:%M:%S.%f") is True
    assert time.is_time_in_given_format('2017-02-14 08:08:12.999', "%Y-%m-%dT%H:%M:%S.%f") is False
