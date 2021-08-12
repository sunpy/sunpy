import pytest

from astropy.time import Time, TimeUnixTai
from astropy.time.formats import erfa

# This registers the TimeTAICDS format with astropy
from sunpy.time import TimeUnixTai58  # NOQA


def test_time_t0():
    """
    Test that the offset at the epoch is zero. Note that the offset
    will be zero in both TAI and UTC
    """
    t = Time('1958-01-01 00:00:00', format='iso', scale='tai')
    assert t.tai.unix_tai_58 == 0.0
    with pytest.warns(erfa.ErfaWarning, match='dubious year'):
        assert t.utc.unix_tai_58 == 0.0


def test_tai_utc_offset():
    """
    Test that the offset between TAI and UTC is 10 s in 1972-01-01, when UTC
    was defined
    """
    t1 = Time('1972-01-01 00:00:00', scale='tai', format='iso')
    t2 = Time('1972-01-01 00:00:00', scale='utc', format='iso')
    assert t2.unix_tai_58 - t1.unix_tai_58 == 10.0


@pytest.mark.parametrize('time', [
    Time('1958-01-01T00:00:00', format='isot', scale='tai'),
    Time('1972-01-01T00:00:00', format='isot', scale='tai'),
    Time('2015-10-25T05:24:08', format='isot', scale='tai'),
    Time('2018-09-17T19:46:25', format='isot', scale='tai'),
])
def test_roundtrip(time):
    time_from_tai_seconds = Time(time.unix_tai_58, scale='tai', format='unix_tai_58')
    assert time_from_tai_seconds == time


@pytest.mark.parametrize('time', [
    Time('1958-01-01T00:00:00', format='isot', scale='tai'),
    Time('1972-01-01T00:00:00', format='isot', scale='tai'),
    Time('2015-10-25T05:24:08', format='isot', scale='tai'),
    Time('2018-09-17T19:46:25', format='isot', scale='tai'),
])
def test_unix_tai_consistency(time):
    diff_unix_tai = Time(TimeUnixTai.epoch_val, scale='tai', format=TimeUnixTai.epoch_format,).unix_tai_58
    assert time.unix_tai + diff_unix_tai == time.unix_tai_58
