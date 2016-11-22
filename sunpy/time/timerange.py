from __future__ import absolute_import, division, print_function

from datetime import timedelta
from datetime import datetime

import astropy.units as u

from sunpy.time import parse_time
from sunpy import config
from sunpy.extern.six.moves import range

TIME_FORMAT = config.get('general', 'time_format')

__all__ = ['TimeRange']


class TimeRange(object):
    """
    An object to handle time ranges.

    .. note::

       Regardless of how a TimeRange is constructed it will always provide a
       positive time range where the start time is before the end time.

    Parameters
    ----------
    a : str, number, `datetime.datetime`
        A time (usually the start time) specified as a parse_time-compatible
        time string, number, or a datetime object.
    b : str, number, `datetime.datetime`, `datetime.timedelta`, `astropy.units.Quantity` (time)
        Another time (usually the end time) specified as a
        parse_time-compatible time string, number, or a datetime object.
        May also be the size of the time range specified as a timedelta object,
        or a `astropy.units.Quantity`.

    Examples
    --------
    >>> from sunpy.time import TimeRange
    >>> time_range = TimeRange('2010/03/04 00:10', '2010/03/04 00:20')
    >>> time_range = TimeRange(('2010/03/04 00:10', '2010/03/04 00:20'))
    >>> import astropy.units as u
    >>> time_range = TimeRange('2010/03/04 00:10', 400 * u.s)
    >>> time_range = TimeRange('2010/03/04 00:10', 400 * u.day)
    """
    def __init__(self, a, b=None):
        """Creates a new TimeRange instance"""
        # If a is a TimeRange object, copy attributes to new instance.
        self._t1 = None
        self._t2 = None

        if isinstance(a, TimeRange):
            self.__dict__ = a.__dict__.copy()
            return

        # Normalize different input types
        if b is None:
            x = parse_time(a[0])
            if len(a) != 2:
                raise ValueError('If b is None a must have two elements')
            else:
                y = a[1]
        else:
            x = parse_time(a)
            y = b

        if isinstance(y, u.Quantity):
            y = timedelta(seconds=y.to('s').value)

        # Timedelta
        if isinstance(y, timedelta):
            if y.days >= 0:
                self._t1 = x
                self._t2 = x + y
            else:
                self._t1 = x + y
                self._t2 = x
            return

        # Otherwise, assume that the second argument is parse_time-compatible
        y = parse_time(y)

        if isinstance(y, datetime):
            if x < y:
                self._t1 = x
                self._t2 = y
            else:
                self._t1 = y
                self._t2 = x

    @property
    def start(self):
        """
        Get the start time

        Returns
        -------
        start : `datetime.datetime`
        """
        return self._t1

    @property
    def end(self):
        """
        Get the end time

        Returns
        -------
        end : `datetime.datetime`
        """
        return self._t2

    @property
    def dt(self):
        """
        Get the length of the time range. Always a positive value.

        Returns
        -------
        dt : `datetime.timedelta`
        """
        return self._t2 - self._t1

    @property
    def center(self):
        """
        Gets the center of the TimeRange instance.

        Returns
        -------
        value : `datetime.datetime`
        """
        return self.start + self.dt // 2

    @property
    def hours(self):
        """
        Get the number of hours elapsed.

        Returns
        -------
        value : `astropy.units.Quantity`
        """
        return self._duration.to('hour')

    @property
    def days(self):
        """
        Gets the number of days elapsed.

        Returns
        -------
        value : `astropy.units.Quantity`
        """
        return self._duration.to('d')

    @property
    def seconds(self):
        """
        Gets the number of seconds elapsed.

        Returns
        -------
        value : `astropy.units.Quantity`
        """
        return self._duration.to('s')

    @property
    def minutes(self):
        """
        Gets the number of minutes elapsed.

        Returns
        -------
        value : `astropy.units.Quantity`
        """
        return self._duration.to('min')

    @property
    def _duration(self):
        """
        The duration of the time range.

        Returns
        -------
        value : `astropy.units.Quantity`
        """
        result = self.dt.microseconds * u.Unit('us') + self.dt.seconds * u.Unit('s') + self.dt.days * u.Unit('day')
        return result

    def __repr__(self):
        """
        Returns a human-readable representation of the TimeRange instance."""

        t1 = self.start.strftime(TIME_FORMAT)
        t2 = self.end.strftime(TIME_FORMAT)
        center = self.center.strftime(TIME_FORMAT)

        return ('    Start:'.ljust(11) + t1 +
                '\n    End:'.ljust(12) + t2 +
                '\n    Center:'.ljust(12) + center +
                '\n    Duration:'.ljust(12) + str(self.days.value) + ' days or' +
                '\n    '.ljust(12) + str(self.hours.value) + ' hours or' +
                '\n    '.ljust(12) + str(self.minutes.value) + ' minutes or' +
                '\n    '.ljust(12) + str(self.seconds.value) + ' seconds' +
                '\n')

    def split(self, n=2):
        """
        Splits the TimeRange into multiple equally sized parts.

        Parameters
        ----------
        n : int
            The number of times to split the time range (must > 1)

        Returns
        -------
        time ranges: list
            An list of equally sized TimeRange objects between
            the start and end times.

        Raises
        ------
        ValueError
            If requested amount is less than 1
        """
        if n <= 0:
            raise ValueError('n must be greater than or equal to 1')
        subsections = []
        previous_time = self.start
        next_time = None
        for _ in range(n):
            next_time = previous_time + self.dt // n
            next_range = TimeRange(previous_time, next_time)
            subsections.append(next_range)
            previous_time = next_time
        return subsections

    def window(self, cadence, window):
        """
        Split the TimeRange up into a series of TimeRange windows,
        'window' long, between the start and end with a cadence of 'cadence'.

        Parameters
        ----------
        cadence : `astropy.units.Quantity`, `datetime.timedelta`
            Cadence in seconds or a timedelta instance
        window : `astropy.units.quantity`, `datetime.timedelta`
            The length of the Time's, assumed to be seconds if int.

        Returns
        -------
        time ranges : list
            A list of TimeRange objects, that are window long and separated by
            cadence.

        Examples
        --------
        >>> import astropy.units as u
        >>> from sunpy.time import TimeRange
        >>> time_range = TimeRange('2010/03/04 00:10', '2010/03/04 01:20')
        >>> time_range.window(60*60*u.s, window=12*u.s)   # doctest: +NORMALIZE_WHITESPACE
        [    Start: 2010-03-04 00:10:00
            End:   2010-03-04 00:10:12
            Center:2010-03-04 00:10:06
            Duration:0.000138888888889 days or
                   0.00333333333333 hours or
                   0.2 minutes or
                   12.0 seconds
         ,    Start: 2010-03-04 01:10:00
            End:   2010-03-04 01:10:12
            Center:2010-03-04 01:10:06
            Duration:0.000138888888889 days or
                   0.00333333333333 hours or
                   0.2 minutes or
                   12.0 seconds
         ,    Start: 2010-03-04 02:10:00
            End:   2010-03-04 02:10:12
            Center:2010-03-04 02:10:06
            Duration:0.000138888888889 days or
                   0.00333333333333 hours or
                   0.2 minutes or
                   12.0 seconds
          ]
        """
        if not isinstance(window, timedelta):
            window = timedelta(seconds=window.to('s').value)
        if not isinstance(cadence, timedelta):
            cadence = timedelta(seconds=cadence.to('s').value)

        n = 1
        times = [TimeRange(self.start, self.start + window)]
        while times[-1].end < self.end:
            times.append(TimeRange(self.start + cadence * n,
                                   self.start + cadence * n + window))
            n += 1
        return times

    def next(self):
        """Shift the time range forward by the amount of time elapsed"""
        dt = self.dt
        self._t1 = self._t1 + dt
        self._t2 = self._t2 + dt

        return self

    def previous(self):
        """Shift the time range backward by the amount of time elapsed"""
        dt = self.dt
        self._t1 = self._t1 - dt
        self._t2 = self._t2 - dt

        return self

    def extend(self, dt_start, dt_end):
        """Extend the time range forwards and backwards

        Parameters
        ----------
        dt_start : `datetime.timedelta`
            The amount to shift the start time
        dt_end : `datetime.timedelta`
            The amount to shift the end time
        """
        # Only a timedelta object is acceptable here
        self._t1 = self._t1 + dt_start
        self._t2 = self._t2 + dt_end

    def __contains__(self, time):
        """
        Checks whether the given time lies within this range.
        Both limits are inclusive (i.e. __contains__(t1) and __contains__(t2)
        always return true)

        Parameters
        ----------
        time : `datetime.datetime`, str
            A parse_time-compatible time to be checked.

        Returns
        -------
        value : bool
            True if time lies between start and end, False otherwise.

        Examples
        --------
        >>> from sunpy.time import TimeRange
        >>> time1 = '2014/5/5 12:11'
        >>> time2 = '2012/5/5 12:11'
        >>> time_range = TimeRange('2014/05/04 13:54', '2018/02/03 12:12')
        >>> time1 in time_range
        True
        >>> time2 in time_range
        False
        """
        this_time = parse_time(time)
        return this_time >= self.start and this_time <= self.end
