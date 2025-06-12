import pytest

from asdf.testing.helpers import roundtrip_object
from astropy.time import Time

# This registers TimeUTime in astropy.time.Time
from sunpy.time import TimeUTime  # NOQA

pytest.importorskip("asdf_astropy", "0.8.0")


def assert_roundtrip_time(old):
    new = roundtrip_object(old)
    assert new.utime == old.utime


def test_utime():
    assert_roundtrip_time(Time(1231853576.0, format='utime'))


def test_tai_seconds():
    assert_roundtrip_time(Time(1231853576.0, format='tai_seconds'))
