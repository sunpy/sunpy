from __future__ import absolute_import

from datetime import timedelta
from datetime import datetime

from sunpy.time import parse_time

__all__ = ['TimeRange']

class TimeRange:
    """
    Timerange(a, b) or Timerange((a, b))

    An object to handle time ranges.

    Parameters
    ----------
    a : the start time specified as a time string, or datetime object
        A 2d list or ndarray containing the map data
    b : the end time specified as a time string or datetime object
        or the length of the time range specified as a timedelta object, or
        number of seconds

    Attributes
    ----------
    t1 : datetime
        The start time of the time range
    t2 : datetime
        The end time of the time range
    dt : timediff
        The difference in time between the start time and end time

    Examples
    --------
    >>> time_range = TimeRange('2010/03/04 00:10', '2010/03/04 00:20')

    >>> time_range = TimeRange('2010/03/04 00:10', 400)

    References
    ----------
    | http://docs.scipy.org/doc/numpy/reference/arrays.classes.html

    """
    def __init__(self, a, b=None):
        """Creates a new TimeRange instance"""
        # If a is a TimeRange object, copy attributes to new instance.
        if isinstance(a, TimeRange):
            self.__dict__ = a.__dict__.copy()
            return

        # Normalize different input types
        if b is None:
            x = a[0]
            y = a[1]
        else:
            x = a
            y = b

        # Start time
        self.t1 = parse_time(x)

        # End date
        if isinstance(y, str):
            self.t2 = parse_time(y)

        # Datetime
        if isinstance(y,datetime):
            self.t2 = y

        # Timedelta
        if isinstance(y, timedelta):
            self.t2 = self.t1 + y

        # Seconds offset
        if isinstance(y, (float, int)):
            self.t2 = self.t1 + timedelta(0, y)

        self.dt = self.t2 - self.t1

    def __repr__(self):
        """Returns a human-readable representation of the TimeRange instance."""
        TIME_FORMAT = "%Y/%m/%d %H:%M:%S"

        t1 = self.t1.strftime(TIME_FORMAT)
        t2 = self.t2.strftime(TIME_FORMAT)
        center = self.center().strftime(TIME_FORMAT)

        return ('\tStart:'.ljust(11) + t1 +
                '\n\tEnd:'.ljust(12) + t2 +
                '\n\tCenter:'.ljust(12) + center +
                '\n\tDuration:'.ljust(12) + str(self.days()) + ' days or' +
                '\n\t'.ljust(12) +  str(self.minutes()) + ' minutes or' +
                '\n\t'.ljust(12) +  str(self.seconds()) + ' seconds' +
                '\n')

    def center(self):
        """Gets the center of the TimeRange instance"""
        return self.t1 + self.dt / 2

    def split(self,n=2):
        """Splits the TimeRange into multiple equally sized parts

        Accepts a value greater than or equal to 1 as input, and
        returns an array of equally sized TimeRange objects between
        t1 and t2.

        Raises a ValueError if requested amount is less than 1

        """

        if n <= 0:
            raise ValueError('n must be greater than or equal to 1')
        subsections = []
        previous_time = self.start()
        next_time = None
        for _ in range(n):
            next_time = previous_time + self.dt/n
            next_range = TimeRange(previous_time,next_time)
            subsections.append(next_range)
            previous_time = next_time
        return subsections
    
    def days(self):
        """Gets the number of days elapsed."""
        return self.dt.days

    def start(self):
        """Gets the start date"""
        return self.t1

    def end(self):
        """Gets the start date"""
        return self.t2

    def seconds(self):
        """Gets the number of seconds elapsed."""
        return (self.dt.microseconds + (self.dt.seconds + self.dt.days * 24 * 3600) * 1e6) / 1e6

    def minutes(self):
        """Gets the number of minutes elapsed."""
        return self.seconds() / 60.0

    def next(self):
        """Shift the time range forward by the amount of time elapsed"""
        self.t1 = self.t1 + self.dt
        self.t2 = self.t2 + self.dt

        return self

    def previous(self):
        """Shift the time range backward by the amount of time elapsed"""
        self.t1 = self.t1 - self.dt
        self.t2 = self.t2 - self.dt

        return self

    def extend(self, t_backwards, t_forwards):
        """Extend the time range forwards and backwards by arbitrary amounts"""
        # Only a timedelta object is acceptable here
        self.t1 = self.t1 + t_backwards
        self.t2 = self.t2 + t_forwards