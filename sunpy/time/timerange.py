"""
This module provies a object that can handle a time range.
"""
from datetime import timedelta

import astropy.units as u
from astropy.time import Time, TimeDelta

from sunpy import config
from sunpy.time import is_time_equal, parse_time
from sunpy.time.time import _variables_for_parse_time_docstring
from sunpy.util.decorators import add_common_docstring

TIME_FORMAT = config.get('general', 'time_format')

__all__ = ['TimeRange']


@add_common_docstring(**_variables_for_parse_time_docstring())
class TimeRange:
    """
    A class to create and handle time ranges.

    .. note::

       Regardless of how a `sunpy.time.TimeRange` is constructed it will always
       provide a positive time range where the start time is before the end time.

       ``__contains__`` has been implemented which means you can
       check if a time is within the time range you have created.
       Please see the example section below.

    Parameters
    ----------
    a : {parse_time_types}
        A time (the start time) specified as a parse_time-compatible
        time string, number, or a datetime object.
    b : {parse_time_types}
        Another time (the end time) specified as a
        parse_time-compatible time string, number, or a datetime object.
        May also be the size of the time range specified as a timedelta object,
        or a `~astropy.units.Quantity`.

    Examples
    --------
    >>> from sunpy.time import TimeRange
    >>> time_range = TimeRange('2010/03/04 00:10', '2010/03/04 00:20')
    >>> time_range = TimeRange(('2010/03/04 00:10', '2010/03/04 00:20'))
    >>> import astropy.units as u
    >>> time_range = TimeRange('2010/03/04 00:10', 400 * u.s)
    >>> TimeRange('2010/03/04 00:10', 400 * u.day)
       <sunpy.time.timerange.TimeRange object at ...>
        Start: 2010-03-04 00:10:00
        End:   2011-04-08 00:10:00
        Center:2010-09-20 00:10:00
        Duration:400.0 days or
               9600.0 hours or
               576000.0 minutes or
               34560000.0 seconds
    <BLANKLINE>
    >>> time1 = '2014/5/5 12:11'
    >>> time2 = '2012/5/5 12:11'
    >>> time_range = TimeRange('2014/05/04 13:54', '2018/02/03 12:12')
    >>> time1 in time_range
    True
    >>> time2 in time_range
    False
    """

    def __init__(self, a, b=None, format=None):
        # If a is a TimeRange object, copy attributes to new instance.
        self._t1 = None
        self._t2 = None

        if isinstance(a, TimeRange):
            self.__dict__ = a.__dict__.copy()
            return

        # Normalize different input types
        if b is None:
            x = parse_time(a[0], format=format)
            if len(a) != 2:
                raise ValueError('If b is None a must have two elements')
            else:
                y = a[1]
        else:
            x = parse_time(a, format=format)
            y = b

        if isinstance(y, u.Quantity):
            y = TimeDelta(y)

        if isinstance(y, timedelta):
            y = TimeDelta(y, format='datetime')

        # Timedelta
        if isinstance(y, TimeDelta):
            if y.jd >= 0:
                self._t1 = x
                self._t2 = x + y
            else:
                self._t1 = x + y
                self._t2 = x
            return

        # Otherwise, assume that the second argument is parse_time-compatible
        y = parse_time(y, format=format)

        if isinstance(y, Time):
            if x < y:
                self._t1 = x
                self._t2 = y
            else:
                self._t1 = y
                self._t2 = x

    @property
    def start(self):
        """
        Get the start time.

        Returns
        -------
        `astropy.time.Time`
            The start time.
        """
        return self._t1

    @property
    def end(self):
        """
        Get the end time.

        Returns
        -------
        `astropy.time.Time`
            The end time.
        """
        return self._t2

    @property
    def dt(self):
        """
        Get the length of the time range. Always a positive value.

        Returns
        -------
        `astropy.time.TimeDelta`
            The difference between the start and the end time.
        """
        return self._t2 - self._t1

    @property
    def center(self):
        """
        Gets the center of the time range.

        Returns
        -------
        `astropy.time.Time`
           The center time.
        """
        return self.start + self.dt / 2

    @property
    def hours(self):
        """
        Get the number of hours elapsed.

        Returns
        -------
        `astropy.units.Quantity`
           The amount of hours between the start and end time.
        """
        return self.dt.to('hour')

    @property
    def days(self):
        """
        Gets the number of days elapsed.

        Returns
        -------
        `astropy.units.Quantity`
            The amount of days between the start and end time.
        """
        return self.dt.to('d')

    @property
    def seconds(self):
        """
        Gets the number of seconds elapsed.

        Returns
        -------
        `astropy.units.Quantity`
           The amount of seconds between the start and end time.
        """
        return self.dt.to('s')

    @property
    def minutes(self):
        """
        Gets the number of minutes elapsed.

        Returns
        -------
        `astropy.units.Quantity`
           The amount of minutes between the start and end time.
        """
        return self.dt.to('min')

    def __eq__(self, other):
        """
        Check that two `sunpy.time.TimeRange` have the same start and end
        datetime.

        Parameters
        ----------
        other : `~sunpy.time.timerange.TimeRange`
            The second `sunpy.time.TimeRange` to compare to.

        Returns
        -------
        `bool`
            `True` if equal, `False` otherwise.
        """
        if isinstance(other, TimeRange):
            return is_time_equal(
                self.start, other.start) and is_time_equal(self.end, other.end)

        return NotImplemented

    def __ne__(self, other):
        """
        Check two `sunpy.time.TimeRange` have different start or end datetimes.

        Parameters
        ----------
        other : `~sunpy.time.timerange.TimeRange`
            The second `sunpy.time.TimeRange` to compare to.

        Returns
        -------
        `bool`
            `True` if non-equal, `False` otherwise.
        """
        if isinstance(other, TimeRange):
            return not (is_time_equal(
                self.start, other.start) and is_time_equal(self.end, other.end))

        return NotImplemented

    def __repr__(self):
        """
        Returns a human-readable representation of `sunpy.time.TimeRange`.
        """

        t1 = self.start.strftime(TIME_FORMAT)
        t2 = self.end.strftime(TIME_FORMAT)
        center = self.center.strftime(TIME_FORMAT)
        fully_qualified_name = f'{self.__class__.__module__}.{self.__class__.__name__}'

        return ('   <{} object at {}>'.format(fully_qualified_name, hex(id(self))) +
                '\n    Start:'.ljust(12) + t1 +
                '\n    End:'.ljust(12) + t2 +
                '\n    Center:'.ljust(12) + center +
                '\n    Duration:'.ljust(12) + str(self.days.value) + ' days or' +
                '\n    '.ljust(12) + str(self.hours.value) + ' hours or' +
                '\n    '.ljust(12) + str(self.minutes.value) + ' minutes or' +
                '\n    '.ljust(12) + str(self.seconds.value) + ' seconds' +
                '\n')

    def split(self, n=2):
        """
        Splits the time range into multiple equally sized parts.

        Parameters
        ----------
        n : `int`, optional
            The number of times to split the time range (must >= 1).
            Defaults to 2.

        Returns
        -------
        `list`
            A list of equally sized `sunpy.time.TimeRange` between the start and end times.
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
        Split the time range up into a series of `~sunpy.time.TimeRange` that
        are ``window`` long, with a cadence of ``cadence``.

        Parameters
        ----------
        cadence : `astropy.units.Quantity`, `astropy.time.TimeDelta`
            Cadence.
        window : `astropy.units.quantity`, `astropy.time.TimeDelta`
            The length of window.

        Returns
        -------
        `list`
            A list of `~sunpy.time.TimeRange`, that are ``window`` long
            and separated by ``cadence``.

        Examples
        --------
        >>> import astropy.units as u
        >>> from sunpy.time import TimeRange
        >>> time_range = TimeRange('2010/03/04 00:10', '2010/03/04 01:20')
        >>> time_range.window(60*60*u.s, window=12*u.s)   # doctest:  +SKIP
        [   <sunpy.time.timerange.TimeRange object at 0x7f0214bfc208>
            Start: 2010-03-04 00:10:00
            End:   2010-03-04 00:10:12
            Center:2010-03-04 00:10:06
            Duration:0.0001388888888888889 days or
                    0.003333333333333333 hours or
                    0.2 minutes or
                    12.0 seconds,
            <sunpy.time.timerange.TimeRange object at 0x7f01fe43ac50>
            Start: 2010-03-04 01:10:00
            End:   2010-03-04 01:10:12
            Center:2010-03-04 01:10:06
            Duration:0.0001388888888888889 days or
                    0.003333333333333333 hours or
                    0.2 minutes or
                    12.0 seconds,
            <sunpy.time.timerange.TimeRange object at 0x7f01fb90b898>
            Start: 2010-03-04 02:10:00
            End:   2010-03-04 02:10:12
            Center:2010-03-04 02:10:06
            Duration:0.0001388888888888889 days or
                    0.003333333333333333 hours or
                    0.2 minutes or
                    12.0 seconds]
        """
        # TODO: After astropy 3.1 remove this check
        if isinstance(window, timedelta):
            window = TimeDelta(window, format="datetime")
        if isinstance(cadence, timedelta):
            cadence = TimeDelta(cadence, format="datetime")

        if not isinstance(window, TimeDelta):
            window = TimeDelta(window)
        if not isinstance(cadence, TimeDelta):
            cadence = TimeDelta(cadence)

        n = 1
        times = [TimeRange(self.start, self.start + window)]
        while times[-1].end < self.end:
            times.append(TimeRange(self.start + cadence * n,
                                   self.start + cadence * n + window))
            n += 1
        return times

    def next(self):
        """
        Shift the time range forward by the amount of time elapsed.
        """
        dt = self.dt
        self._t1 = self._t1 + dt
        self._t2 = self._t2 + dt

        return self

    def previous(self):
        """
        Shift the time range backward by the amount of time elapsed.
        """
        dt = self.dt
        self._t1 = self._t1 - dt
        self._t2 = self._t2 - dt

        return self

    def extend(self, dt_start, dt_end):
        """
        Extend the time range forwards and backwards.

        Parameters
        ----------
        dt_start : `astropy.time.TimeDelta`
            The amount to shift the start time.
        dt_end : `astropy.time.TimeDelta`
            The amount to shift the end time.
        """
        # TODO: Support datetime.timedelta
        self._t1 = self._t1 + dt_start
        self._t2 = self._t2 + dt_end

    def get_dates(self):
        """
        Return all partial days contained within the time range.
        """
        delta = self.end.to_datetime().date() - self.start.to_datetime().date()
        dates = []
        dates = [
            parse_time(self.start.strftime('%Y-%m-%d')) + TimeDelta(i*u.day)
            for i in range(delta.days + 1)
        ]
        return dates

    @add_common_docstring(**_variables_for_parse_time_docstring())
    def __contains__(self, time):
        """
        Checks whether the given time lies within this range. Both limits are
        inclusive (i.e., ``__contains__(t1)`` and ``__contains__(t2)`` always
        return `True).that.

        Parameters
        ----------
        time : {parse_time_types}
            {parse_time_desc}

        Returns
        -------
        `bool`
            `True` if time lies between start and end, `False` otherwise.

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
