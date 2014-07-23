from __future__ import absolute_import

from datetime import timedelta
from datetime import datetime

from sunpy.time import parse_time

__all__ = ['TimeRange']


class TimeRange(object):
    """
    An object to handle time ranges.

    Parameters
    ----------
    a : str, number, datetime,
        A time (usually the start time) specified as a parse_time-compatible
        time string or number, or a datetime object.
    b : str, number, timedelta
        Another time (usually the end time) specified as a
        parse_time-compatible time string, or a datetime object.
        May also be the size of the time range specified as a timedelta object,
        or the number of seconds (positive or negative)

    Properties
    ----------
    start : datetime
        The start time of the time range (always the smaller time)
    end : datetime
        The end time of the time range (always the larger time)

    Attributes
    ----------
    dt : timediff
        The difference in time between the start and end time. Always a
        positive value.

    Examples
    --------
    >>> from sunpy.time import TimeRange
    >>> time_range = TimeRange('2010/03/04 00:10', '2010/03/04 00:20')
    >>> time_range = TimeRange(('2010/03/04 00:10', '2010/03/04 00:20'))
    >>> time_range = TimeRange('2010/03/04 00:10', 400)
    >>> time_range = TimeRange('2010/03/04 00:10', -400)

    References
    ----------
    | http://docs.scipy.org/doc/numpy/reference/arrays.classes.html

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

        if isinstance(y, str):
            y = parse_time(y)

        if isinstance(y, datetime):
            if x < y:
                self._t1 = x
                self._t2 = y
            else:
                self._t1 = y
                self._t2 = x

        if isinstance(y, (float, int)):
            y = timedelta(seconds=y)

        # Timedelta
        if isinstance(y, timedelta):
            if y.days >= 0:
                self._t1 = x
                self._t2 = x + y
            else:
                self._t1 = x + y
                self._t2 = x

    @property
    def start(self):
        """Get the start time

        Returns
        -------
        start : datetime
        """
        return self._t1

    @property
    def end(self):
        """Get the end time

        Returns
        -------
        end : datetime
        """
        return self._t2

    @property
    def dt(self):
        """Get the length of the time range

        Returns
        -------
        dt : timedelta
        """
        return self._t2 - self._t1

    def __repr__(self):
        """Returns a human-readable representation of the TimeRange
        instance."""
        TIME_FORMAT = "%Y/%m/%d %H:%M:%S"

        t1 = self.start.strftime(TIME_FORMAT)
        t2 = self.end.strftime(TIME_FORMAT)
        center = self.center().strftime(TIME_FORMAT)

        return ('    Start:'.ljust(11) + t1 +
                '\n    End:'.ljust(12) + t2 +
                '\n    Center:'.ljust(12) + center +
                '\n    Duration:'.ljust(12) + str(self.days()) + ' days or' +
                '\n    '.ljust(12) + str(self.minutes()) + ' minutes or' +
                '\n    '.ljust(12) + str(self.seconds()) + ' seconds' +
                '\n')

    def center(self):
        """Gets the center of the TimeRange instance

        Returns
        -------
        value : datetime
        """
        return self.start + self.dt / 2

    def split(self, n=2):
        """Splits the TimeRange into multiple equally sized parts

        Parameters
        ----------
        n : int
            The number of times to split the time range
            (must be an integer greater than 1)

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
            next_time = previous_time + self.dt / n
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
        cadence: int or timedelta
            Cadence in seconds or a timedelta instance
        window: int or timedelta
            The length of the Time's, assumed to be seconds if int.

        Returns
        -------
        time ranges: list
            A list of TimeRange objects, that are window long and seperated by
            cadence.

        Examples
        --------
        # To get one 12 second long window every hour within the timerange:
        >>> TimeRange.window(60*60, window=12)
        """
        if not isinstance(window, timedelta):
            window = timedelta(seconds=window)
        if not isinstance(cadence, timedelta):
            cadence = timedelta(seconds=cadence)

        n = 1
        times = [TimeRange(self.start, self.start + window)]
        while times[-1].end < self.end:
            times.append(TimeRange(self.start + cadence * n,
                                   self.start + cadence * n + window))
            n += 1
        return times

    def days(self):
        """Gets the number of days elapsed."""
        return self.dt.days

    def seconds(self):
        """Gets the number of seconds elapsed."""
        return (self.dt.microseconds +
               (self.dt.seconds + self.dt.days * 24 * 3600) * 1e6) / 1e6

    def minutes(self):
        """Gets the number of minutes elapsed."""
        return self.seconds() / 60.0

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
        dt_start : timedelta
            The amount to shift the start time
        dt_end : timedelta
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
        time: datetime or str
            A parse_time-compatible time to be checked.

        Returns
        -------
        value : bool
            True if time lies between start and end, False otherwise.

        Example
        -------
        >>> time_range = TimeRange('2014/05/04 13:54', '2018/02/03 12:12')
        >>> time in time_range
        """
        t = parse_time(time)
        return t >= self.start and t <= self.end
