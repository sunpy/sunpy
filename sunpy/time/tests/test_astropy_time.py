from sunpy.time.astropy_time import Time

import numpy as np

import pytest


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
