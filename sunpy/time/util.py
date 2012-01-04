from __future__ import absolute_import

__all__ = ["TimeRange"]

import datetime
from sunpy.util import anytim

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
        The end time of the time range
    center: datetime
        The center of the time range
    dt : timediff
        The difference in time between the start time and end time
    show : str
        Display the time string in a human readable format
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
            if isinstance(a[0],str) or isinstance(a[0], float):
                self.t1 = anytim(a[0])
            if isinstance(a[1],str) or isinstance(a[1], float):
                self.t2 = anytim(a[1])
            if type(a[1]) == type(datetime.timedelta(1)):
                self.t2 = self.t1 + a[1]
            if isinstance(a[1], int):
                self.t2 = self.t1 + datetime.timedelta(0,a[1])                
        else:            
            if isinstance(a,str) or isinstance(a, float):
                self.t1 = anytim(a)
            if isinstance(b,str) or isinstance(b, float):
                self.t2 = anytim(b)
            if type(b) == type(datetime.timedelta(1)):
                self.t2 = self.t1 + b
            if isinstance(b, int):
                self.t2 = self.t1 + datetime.timedelta(0,b) 
        self.dt = self.t2 - self.t1
    
    def center(self):
        return self.t1 + self.dt/2
    
    def days(self):
        return self.dt.days
    
    def seconds(self):
        return self.dt.total_seconds()
    
    def minutes(self):
        return self.dt.total_seconds()/60.0
    
    def show(self):
        print(self.t1.strftime("%Y/%m/%d %H:%M:%S") + ' ' + 
              self.t2.strftime("%Y/%m/%d %H:%M:%S"))