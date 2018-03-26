from __future__ import absolute_import

from astropy.time.formats import erfa, TimeFromEpoch

__all__ = ['TimeUTime']


class TimeUTime(TimeFromEpoch):
    """
    Seconds from 1979-01-01 00:00:00 UTC. Same as Unix time
    but this starts 9 years later.

    This time format is included for historical reasons. Some
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
