from __future__ import absolute_import
__all__ = ["TimeRange"]

class TimeRange:
    """
    Timerange(a, b) or Timerange((a,b))

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
    days : float
        Number of days in the time range
    minutes: float
        Number of minutes in the time range
    seconds: float
        Number of seconds in the time range
    next : None
        Shift the start time (t1) and end time (t2) by adding dt in 
        the current instance
    previous : None
        Shift the start time (t1) and end time (t2) by subtracting dt in 
        the current instance
   
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
        from sunpy.time import parse_time

        if b is None:
            if isinstance(a[0],str) or isinstance(a[0], float):
                self.t1 = parse_time(a[0])
            if isinstance(a[1],str) or isinstance(a[1], float):
                self.t2 = parse_time(a[1])
            if type(a[1]) == type(timedelta(1)):
                self.t2 = self.t1 + a[1]
            if isinstance(a[1], int):
                self.t2 = self.t1 + timedelta(0,a[1])                
        else:            
            if isinstance(a,str) or isinstance(a, float):
                self.t1 = parse_time(a)
            if isinstance(b,str) or isinstance(b, float):
                self.t2 = parse_time(b)
            if type(b) == type(timedelta(1)):
                self.t2 = self.t1 + b
            if isinstance(b, int):
                self.t2 = self.t1 + timedelta(0,b) 
        self.dt = self.t2 - self.t1
    
    def __repr__(self):
        TIME_FORMAT = "%Y/%m/%d %H:%M:%S"
        print('\tStart time: ' + self.t1.strftime(TIME_FORMAT) + '\n' + 
              '\tEnd time: ' + self.t2.strftime(TIME_FORMAT) + '\n' + 
              '\tCenter time: ' + self.center().strftime(TIME_FORMAT) + '\n' + 
              '\tDuration: ' + str(self.days()) + ' days or\n\t' + 
                            str(self.minutes()) + ' min or\n\t' + 
                            str(self.seconds()) + ' seconds')

    def center(self):
        return self.t1 + self.dt/2
    
    def days(self):
        return self.dt.days
    
    def seconds(self):
        return self.dt.total_seconds()
    
    def minutes(self):
        return self.dt.total_seconds()/60.0
    
    def next(self):
        self.t1 = self.t1 + self.dt
        self.t2 = self.t2 + self.dt
    
    def previous(self):
        self.t1 = self.t1 - self.dt
        self.t2 = self.t2 - self.dt