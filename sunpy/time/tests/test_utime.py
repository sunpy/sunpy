from astropy.time import Time

# This registers TimeUTime in astropy.time.Time
from sunpy.time import TimeUTime  # NOQA


def test_utime_t0():
    assert Time('1979-01-01T00:00:00').utime == 0.0


def test_utime_random_date():
    assert Time('2018-01-13T13:32:56').utime == 1231853576.0


def test_conversion_from_utime():
    # This is expected to fail until astropy issue 7092 is resolved.
    t2 = Time(1231853576.0, format='utime')
    assert t2.isot == '2018-01-13T13:32:56.000'
