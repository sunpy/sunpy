from datetime import datetime, timedelta

import numpy as np
import pytest

from astropy.time import TimeDelta
from astropy.time import Time
import astropy.units as u

from sunpy.time import is_time_equal


def test_strftime_scalar():
    """Test of Time.strftime
    """
    time_string = '2010-09-03 06:00:00'
    t = Time(time_string)

    for format in t.FORMATS:
        t.format = format
        assert t.strftime('%Y-%m-%d %H:%M:%S') == time_string


def test_strftime_array():
    tstrings = ['2010-09-03 00:00:00', '2005-09-03 06:00:00', '1995-12-31 23:59:60']
    t = Time(tstrings)

    for format in t.FORMATS:
        t.format = format
        assert t.strftime('%Y-%m-%d %H:%M:%S').tolist() == tstrings


def test_strftime_array_2():
    tstrings = [['1998-01-01 00:00:01', '1998-01-01 00:00:02'],
                ['1998-01-01 00:00:03', '1995-12-31 23:59:60']]
    tstrings = np.array(tstrings)

    t = Time(tstrings)

    for format in t.FORMATS:
        t.format = format
        assert np.all(t.strftime('%Y-%m-%d %H:%M:%S') == tstrings)
        assert t.strftime('%Y-%m-%d %H:%M:%S').shape == tstrings.shape


def test_strftime_leapsecond():
    time_string = '1995-12-31 23:59:60'
    t = Time(time_string)

    for format in t.FORMATS:
        t.format = format
        assert t.strftime('%Y-%m-%d %H:%M:%S') == time_string


def test_strptime_scalar():
    """Test of Time.strptime
    """
    time_string = '2007-May-04 21:08:12'
    time_object = Time('2007-05-04 21:08:12')
    t = Time.strptime(time_string, '%Y-%b-%d %H:%M:%S')

    assert t == time_object


def test_strptime_array():
    """Test of Time.strptime
    """
    tstrings = [['1998-Jan-01 00:00:01', '1998-Jan-01 00:00:02'],
                ['1998-Jan-01 00:00:03', '1998-Jan-01 00:00:04']]
    tstrings = np.array(tstrings)

    time_object = Time([['1998-01-01 00:00:01', '1998-01-01 00:00:02'],
                        ['1998-01-01 00:00:03', '1998-01-01 00:00:04']])
    t = Time.strptime(tstrings, '%Y-%b-%d %H:%M:%S')

    assert np.all(t == time_object)
    assert t.shape == tstrings.shape


def test_strptime_badinput():
    tstrings = [1, 2, 3]
    with pytest.raises(TypeError):
        Time.strptime(tstrings, '%S')


def test_strptime_input_bytes_scalar():
    time_string = b'2007-May-04 21:08:12'
    time_object = Time('2007-05-04 21:08:12')
    t = Time.strptime(time_string, '%Y-%b-%d %H:%M:%S')

    assert t == time_object


def test_strptime_input_bytes_array():
    tstrings = [[b'1998-Jan-01 00:00:01', b'1998-Jan-01 00:00:02'],
                [b'1998-Jan-01 00:00:03', b'1998-Jan-01 00:00:04']]
    tstrings = np.array(tstrings)

    time_object = Time([['1998-01-01 00:00:01', '1998-01-01 00:00:02'],
                        ['1998-01-01 00:00:03', '1998-01-01 00:00:04']])
    t = Time.strptime(tstrings, '%Y-%b-%d %H:%M:%S')

    assert np.all(t == time_object)
    assert t.shape == tstrings.shape


def test_strptime_leapsecond():
    time_obj1 = Time('1995-12-31T23:59:60', format='isot')
    time_obj2 = Time.strptime('1995-Dec-31 23:59:60', '%Y-%b-%d %H:%M:%S')

    assert time_obj1 == time_obj2


def testis_time_equal():
    t1 = Time('1995-12-31T23:59:60', format='isot')
    t2 = Time('1995-12-31T23:59:60', format='isot')

    assert is_time_equal(t1, t2)

    t1 = Time('1995-12-31T23:59:59', format='isot')
    t2 = Time(datetime(1995, 12, 31, 23, 59, 59), format='datetime')

    assert is_time_equal(t1, t2)

    t1 = Time('1995-12-31T23:59:60', format='isot')
    t2 = Time('1995-12-31T23:59:60', format='isot') + TimeDelta(0*u.day)

    assert is_time_equal(t1, t2)


def testis_time_equal_not_equal():
    t1 = Time('1995-12-31T23:59:59', format='isot')
    t2 = Time('1995-12-31T23:59:60', format='isot')

    assert not is_time_equal(t1, t2)

    t1 = Time('1995-12-31T23:59:59', format='isot')
    t2 = Time('1995-12-31T23:59:59', format='isot') + TimeDelta(2*u.nanosecond)

    assert not is_time_equal(t1, t2)


def test_python_timedelta_scalar():
    td = timedelta(days=1, seconds=1)
    td1 = TimeDelta(td, format='datetime')

    assert td1.sec == 86401.0

    td2 = TimeDelta(86401.0, format='sec')
    assert td2.datetime == td


def test_python_timedelta_vector():
    td = [[timedelta(days=1), timedelta(days=2)],
          [timedelta(days=3), timedelta(days=4)]]

    td1 = TimeDelta(td, format='datetime')

    assert np.all(td1.jd == [[1, 2], [3, 4]])

    td2 = TimeDelta([[1, 2], [3, 4]], format='jd')
    assert np.all(td2.datetime == td)


def test_timedelta_to_datetime():
    td = TimeDelta(1, format='jd')

    assert td.to_datetime() == timedelta(days=1)

    td2 = TimeDelta([[1, 2], [3, 4]], format='jd')
    td = [[timedelta(days=1), timedelta(days=2)],
          [timedelta(days=3), timedelta(days=4)]]

    assert np.all(td2.to_datetime() == td)
