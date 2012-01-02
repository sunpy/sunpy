from __future__ import absolute_import

__all__ = ["TimeRange"]

import datetime
from sunpy.util.util import anytim

class TimeRange:
    """
    Timerange((a, b))

    An object to handle time ranges.
    
    Parameters
    ----------
    a : the start time specified as a time string, or datetime object
        A 2d list or ndarray containing the map data
    b : the end time specified as a time string or datetime object
        or the length of the time range specified as a timediff object, or 
        number of seconds

    Attributes
    ----------
    t1 : datetime
        The start time of the time range
    t2 : datetime
        The end time fo the time range
    dt : timediff
        The difference in time between the start time and end time
    days : float
        Number of days in the time range
    minutes: float
        Number of minutes in the time range
    seconds: float
        Number of seconds in the time range
   
    Examples
    --------
    >>> time_range = TimeRange('2010/03/04 00:10', '2010/03/04 00:20')
    >>> time_range = TimeRange('2010/03/04 00:10', 400)
    >>> time_range = TimeRange(['2010/03/04 00:10', '2010/03/04 00:20'])
    >>> time_range = TimeRange(['2010/03/04 00:10', 400])
    
    See Also
    --------
    
    References
    ----------
    | http://docs.scipy.org/doc/numpy/reference/arrays.classes.html

    """
    def __init__(self, a, b = None):
        if b is None:
            if type(a[0]) == type('string'):
                self.t1 = anytim(a[0])
            if type(a[1]) == type('string'):
                self.t2 = anytim(a[1])
            if type(a[1]) == type(datetime.timedelta(1)):
                self.t2 = self.t1 + a[1]
            if type(a[1]) == type(1):
                self.t2 = self.t1 + datetime.timedelta(0,b)
        else:            
            if type(a) == type('string'):
                self.t1 = anytim(a)
            if type(b) == type('string'):
                self.t2 = anytim(b)
            if type(b) == type(datetime.timedelta(1)):
                self.t2 = self.t1 + b
            if type(b) == type(1):
                self.t2 = self.t1 + datetime.timedelta(0,b)
            
        self.dt = self.t2 - self.t1
            
    def days(self):
        return self.dt.days
    
    def seconds(self):
        return self.dt.total_seconds()
    
    def minutes(self):
        return self.dt.total_seconds()/60.0