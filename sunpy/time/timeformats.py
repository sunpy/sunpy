"""
This module provides time formats specific to solar physics
"""
from astropy.time.formats import TimeFromEpoch, erfa

__all__ = ['TimeUTime', 'TimeUnixTai58']


class TimeUTime(TimeFromEpoch):
    """
    Seconds from 1979-01-01 00:00:00 UTC.

    This is equivalent to `~astropy.time.TimeUnix`, except that the epoch
    is 9 years later. This format is included for historical reasons as some
    people in solar physics prefer using this epoch.

    Examples
    --------
    >>> from astropy.time import Time
    >>> t = Time('2000-01-01T13:53:23')
    >>> print(t.utime)
    662738003.0
    >>> t2 = Time('1979-01-01T00:00:00')
    >>> print(t2.utime)
    0.0
    """
    name = 'utime'
    unit = 1.0 / erfa.DAYSEC  # in days (1 day == 86400 seconds)
    epoch_val = '1979-01-01 00:00:00'
    epoch_val2 = None
    epoch_scale = 'utc'  # Scale for epoch_val class attribute
    epoch_format = 'iso'  # Format for epoch_val class attribute


class TimeUnixTai58(TimeFromEpoch):
    """
    Seconds from 1958-01-01 00:00:00, including leap seconds.

    This is equivalent to `~astropy.time.TimeUnixTai`, except that the
    epoch is 12 years earlier.

    .. note:: For dates on and after 1972-01-01 00:00:00 UTC, this format is
              equivalent to that returned by the `anytim2tai` routine in SSW.
              For dates between 1958-01-01 00:00:00 UTC and 1971-12-31 23:59:59 UTC,
              `anytim2tai` returns a constant difference of 9 s difference between
              UTC and TAI (i.e. the number of leap seconds added since 1958-01-01 00:00:00)
              while `~astropy.time.Time` returns a 0 s  on 1958-01-01 00:00:00 (when
              UT and TAI were synchronized) and increases approximately linearly to a
              differnce of 10 s on 1972-01-01 00:00:00. See the "Known Issues" page
              for a more detailed discussion of this discrepancy.

    Examples
    --------
    >>> from astropy.time import Time
    >>> t = Time('1958-01-01T00:00:00', format='isot', scale='tai')
    >>> t.unix_tai_58
    0.0
    >>> t2 = Time('2015-10-25T05:24:08', format='isot', scale='tai')
    >>> t3 = Time(t2.unix_tai_58, format='unix_tai_58', scale='tai')
    >>> t3
    >>> t3.isot
    '2015-10-25T05:24:08.000'

    References
    ----------
    * `CDS Time Conversion Software README <https://hesperia.gsfc.nasa.gov/ssw/gen/idl/time/aaareadme.txt>`_
    * `anytim2tai routine in SSW <https://hesperia.gsfc.nasa.gov/ssw/gen/idl/time/anytim2tai.pro>`_
    """
    name = 'unix_tai_58'
    unit = 1.0 / erfa.DAYSEC  # in days (1 day == 86400 seconds)
    epoch_val = '1958-01-01 00:00:00'
    epoch_val2 = None
    epoch_scale = 'tai'
    epoch_format = 'iso'
