import pytest

from astropy.time import Time
from astropy.time.formats import erfa

# This registers the TimeTaiSeconds format with astropy
from sunpy.time import *  # NOQA


def test_time_t0():
    """
    Test that the offset at the epoch is zero. Note that the offset
    will be zero in both TAI and UTC
    """
    t = Time('1958-01-01 00:00:00', format='iso', scale='tai')
    assert t.tai.tai_seconds == 0.0
    with pytest.warns(erfa.ErfaWarning, match='dubious year'):
        assert t.utc.tai_seconds == 0.0


def test_tai_utc_offset():
    """
    Test that the offset between TAI and UTC is 10 s in 1972-01-01, when UTC
    was defined
    """
    t1 = Time('1972-01-01 00:00:00', scale='tai', format='iso')
    t2 = Time('1972-01-01 00:00:00', scale='utc', format='iso')
    assert t2.tai_seconds - t1.tai_seconds == 10.0


@pytest.mark.parametrize('time', [
    Time('1958-01-01T00:00:00', format='isot', scale='tai'),
    Time('1972-01-01T00:00:00', format='isot', scale='tai'),
    Time('2015-10-25T05:24:08', format='isot', scale='tai'),
    Time('2018-09-17T19:46:25', format='isot', scale='tai'),
])
def test_roundtrip(time):
    time_from_tai_seconds = Time(time.tai_seconds, scale='tai', format='tai_seconds')
    assert time_from_tai_seconds == time
