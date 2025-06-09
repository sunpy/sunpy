from asdf.testing.helpers import roundtrip_object
from astropy.time import Time

# This registers TimeUTime in astropy.time.Time
from sunpy.time import TimeUTime  # NOQA


def assert_roundtrip_time(old):
    new = roundtrip_object(old)
    assert new.utime == old.utime


def test_utime():
    assert_roundtrip_time(Time(1231853576.0, format='utime'))
