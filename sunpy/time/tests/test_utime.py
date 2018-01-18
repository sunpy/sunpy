from __future__ import absolute_import, division, print_function

from astropy.time import Time


def test_utime_t0():
    assert Time('1979-01-01T00:00:00').utime == 0.0


def test_utime_random_date():
    assert Time('2018-01-13T13:32:56').utime == 1231853576.0
