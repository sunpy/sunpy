"""
This module provides time formats specific to solar physics
"""
from astropy.time.formats import TimeFromEpoch, erfa

__all__ = ['TimeUTime', 'TimeTaiSeconds']


class TimeUTime(TimeFromEpoch):
    """
    UT seconds from 1979-01-01 00:00:00 UTC, ignoring leap seconds.

    Notes
    -----
    This format "ignores" leap seconds by treating each day as spanned by exactly
    86400 seconds, which means that this format's second is not actually uniform in
    duration. On a day without a leap second, this format's second is equal to an
    SI second. On a day with a leap second, this format's second is larger than an
    SI second by 1/86400 of an SI second.

    This format is very similar to the default output format of the ``anytim``
    routine in SSW in that there are exactly 86400 seconds assigned for each day.
    However, ``anytim`` treats the seconds as always equal to an SI second, and thus
    the 86400 seconds span only the first 86400/86401 of the day, and the leap
    second is skipped over. This results in discrepancies of up to a second on days
    with a leap second.

    This format is equivalent to `~astropy.time.TimeUnix`, except that the epoch is
    9 years later.

    References
    ----------
    * `anytim routine in SSW <https://hesperia.gsfc.nasa.gov/ssw/gen/idl/utplot/anytim.pro>`__

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


class TimeTaiSeconds(TimeFromEpoch):
    """
    SI seconds from 1958-01-01 00:00:00, which includes UTC leap seconds.

    1958-01-01 00:00:00 is the defined time when International Atomic Time (TAI)
    and Universal Time (UT) are synchronized. A TAI second has the same length as
    an SI second, but prior to 1972-01-01, a UT second -- then defined to be 1/86400
    of an Earth day -- grew to be longer than than an SI second. 1972-01-01
    00:00:00 UTC is equal to 1972-01-01 00:00:10 TAI. After 1972-01-01, Coordinated
    Universal Time (UTC) is defined with seconds that are the same length as SI
    seconds, but now leap seconds are occasionally added to UTC so that it stays
    synchronized with Earth days.

    Notes
    -----
    This format is equivalent to the output of the SSW ``anytim2tai`` routine, and
    related routines, for times after 1972-01-01. Be aware that the SSW routines
    are not written to provide valid results for times before 1972-01-01.

    This format is equivalent to `~astropy.time.TimeUnixTai`, except that the epoch
    is 12 years earlier.

    References
    ----------
    * `anytim2tai routine in SSW <https://hesperia.gsfc.nasa.gov/ssw/gen/idl/time/anytim2tai.pro>`__

    Examples
    --------
    >>> from astropy.time import Time
    >>> t = Time('1958-01-01T00:00:00', format='isot', scale='tai')
    >>> t.tai_seconds  # doctest: +SKIP
    0.0
    >>> t2 = Time('2015-10-25T05:24:08', format='isot', scale='tai')
    >>> t2.tai_seconds  # doctest: +SKIP
    1824441848.0
    >>> t3 = Time(t2.tai_seconds, format='tai_seconds')  # scale is automatically TAI
    >>> t3.isot
    '2015-10-25T05:24:08.000'
    """
    name = 'tai_seconds'
    unit = 1.0 / erfa.DAYSEC  # in days (1 day == 86400 seconds)
    epoch_val = '1958-01-01 00:00:00'
    epoch_val2 = None
    epoch_scale = 'tai'
    epoch_format = 'iso'
