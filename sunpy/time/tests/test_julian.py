from __future__ import absolute_import, division, print_function

from datetime import datetime

from numpy.testing import assert_almost_equal
import pytest

from sunpy.time import julian

DATETIME_DATE_1 = datetime(1900, 1, 1, 12, 00, 0)
STRING_DATE_1 = '1900/1/1 12:00:00'

DATETIME_DATE_2 = datetime(1984, 11, 23, 14, 21, 5)
STRING_DATE_2 = '1984/11/23 14:21:05'

DATETIME_DATE_3 = datetime(2174, 2, 11, 5, 2, 0)
STRING_DATE_3 = '2174/02/11 05:02:00'

DATETIME_DATE_4 = datetime(814, 1, 28, 23, 59, 59)
STRING_DATE_4 = '0814/01/28 23:59:59'

LANDING = datetime(1966, 2, 3)

def test__all__():
    """should return __all__"""

    assert julian.__all__ == ["julian_day", "julian_centuries"]

def test_julian_day():
    assert julian.julian_day('1900-01-01 12:00') == 2415021.0
    assert julian.julian_day(LANDING) == 2439159.5
    result = julian.julian_day('2000-03-01 15:30:26')
    assert_almost_equal(result, 2451605.1461111, decimal=3)

def test_julian_day1():
    """should return julian day for date 1"""
    expected_day = 2415021.0
    assert julian.julian_day(DATETIME_DATE_1) == expected_day
    assert julian.julian_day(STRING_DATE_1) == expected_day

def test_julian_day2():
    """should return julian day for date 2"""

    expected_day = 2446028.097974537
    assert julian.julian_day(DATETIME_DATE_2) == expected_day
    assert julian.julian_day(STRING_DATE_2) == expected_day

def test_julian_day3():
    """should return julian day for date 3"""

    expected_day = 2515138.7097222223
    assert julian.julian_day(DATETIME_DATE_3) == expected_day
    assert julian.julian_day(STRING_DATE_3) == expected_day

def test_julian_day4():
    """should return julian day for date 4"""

    expected_day = 2018395.499988426
    assert julian.julian_day(DATETIME_DATE_4) == expected_day
    assert julian.julian_day(STRING_DATE_4) == expected_day

def test_julian_day5():
    """should raise value error when passed empty string"""

    pytest.raises(ValueError, julian.julian_day, '')

def test_julian_day6():
    """should raise value error when passed non-date string"""

    pytest.raises(ValueError, julian.julian_day, 'A lovely bunch of coconuts')

def test_julian_centuries1():
    """should return julian century for date 1"""

    expected_century = 2.7378507871321012e-05
    assert julian.julian_centuries(DATETIME_DATE_1) == expected_century
    assert julian.julian_centuries(STRING_DATE_1) == expected_century

def test_julian_centuries2():
    """should return julian century for date 2"""

    expected_century = 0.8489554544705528
    assert julian.julian_centuries(DATETIME_DATE_2) == expected_century
    assert julian.julian_centuries(STRING_DATE_2) == expected_century

def test_julian_centuries3():
    """should return julian century for date 3"""

    expected_century = 2.741100882196367
    assert julian.julian_centuries(DATETIME_DATE_3) == expected_century
    assert julian.julian_centuries(STRING_DATE_3) == expected_century

def test_julian_centuries4():
    """should return julian century for date 4"""

    expected_century = -10.85898699552564
    assert julian.julian_centuries(DATETIME_DATE_4) == expected_century
    assert julian.julian_centuries(STRING_DATE_4) == expected_century


def test_julian_centuries5():
    """should raise value error when passed empty string"""

    pytest.raises(ValueError, julian.julian_centuries, '')

def test_julian_centuries6():
    """should raise value error when passed non-date string"""

    pytest.raises(ValueError, julian.julian_centuries, 'Are you suggesting coconuts migrate?')
