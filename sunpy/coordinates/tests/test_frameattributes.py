# -*- coding: utf-8 -*-

import datetime

import pytest

from astropy.time import Time

from ..frames import Helioprojective
from ..frameattributes import TimeFrameAttributeSunPy

@pytest.fixture
def attr():
    return TimeFrameAttributeSunPy()

def test_now(attr):
    """ We can't actually test the value independantly """
    result, converted = attr.convert_input('now')

    assert isinstance(result, Time)
    assert converted

def test_none(attr):
    """ We can't actually test the value independantly """
    result, converted = attr.convert_input(None)

    assert result is None
    assert not converted

@pytest.mark.parametrize('input', [Time('2012-01-01 00:00:00'), '2012/01/01T00:00:00',
                                   '20120101000000', '2012/01/01 00:00:00'])
def test_convert(attr, input):
    result, converted = attr.convert_input(input)

    output = Time('2012-01-01 00:00:00')

    assert isinstance(result, Time)
    assert result == output


@pytest.mark.parametrize('input', [Time('2012-01-01 00:00:00'), '2012/01/01T00:00:00',
                                   '20120101000000', '2012/01/01 00:00:00'])
def test_on_frame(input):
    hpc1 = Helioprojective(dateobs=input)

    output = Time('2012-01-01 00:00:00')

    assert isinstance(hpc1.dateobs, Time)
    assert hpc1.dateobs == output

def test_non_string():
    output = datetime.datetime.now()

    hpc1 = Helioprojective(dateobs=output)

    assert isinstance(hpc1.dateobs, Time)
    assert hpc1.dateobs == output

def test_on_frame_error():
    with pytest.raises(ValueError):
        hpc1 = Helioprojective(dateobs='ajshdasjdhk')

def test_on_frame_error2():
    with pytest.raises(ValueError):
        hpc1 = Helioprojective(dateobs=17263871263)

def test_array():
    input = Time(['2012-01-01 00:00:00', '2012-01-01 00:00:05'])
    with pytest.raises(ValueError):
        hpc1 = Helioprojective(dateobs=input)
