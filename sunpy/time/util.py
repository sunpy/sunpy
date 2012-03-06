from __future__ import absolute_import

__all__ = ["TimeRange", "julian_day", "julian_centuries", 
           "day_of_year", "break_time"]

# The number of days between Jan 1 1900 and the Julian reference date of 
# 12:00 noon Jan 1, 4713 BC
JULIAN_DAY_ON_NOON01JAN1900 = 2415021.0

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


def julian_day(t=None):
    """Returns the (fractional) Julian day defined as the number of days 
    between the queried day and the reference date of 12:00 (noon) Jan 1, 4713 
    BC."""
    # Good online reference for fractional julian day
    # http://www.stevegs.com/jd_calc/jd_calc.htm
    
    JULIAN_REF_DAY = parse_time('1900/1/1 12:00:00')
    time = parse_time(t)
    
    tdiff = time - JULIAN_REF_DAY
    
    julian = tdiff.days + JULIAN_DAY_ON_NOON01JAN1900
   
    result = julian + 1/24.*(time.hour + time.minute/60.0 + 
                             time.second/(60.*60.))

    # This is because the days in datetime objects start at 00:00, 
    # not 12:00 as for Julian days.
    if time.hour >= 12:
        result = result - 0.5
    else:
        result = result + 0.5

    return result

def julian_centuries(t=None):
    """Returns the number of Julian centuries since 1900 January 0.5."""
    DAYS_IN_YEAR = 36525.0

    result = (julian_day(t) - JULIAN_DAY_ON_NOON01JAN1900) / DAYS_IN_YEAR
    return result

def day_of_year(t=None):
    """Returns the day of year."""
    SECONDS_IN_DAY = 60*60*24.0
    time = parse_time(t)
    time_diff = parse_time(t) - datetime(time.year, 1, 1, 0, 0, 0)
    result = time_diff.days + time_diff.seconds/SECONDS_IN_DAY
    return result

def break_time(t=None):
    """Given a time returns a string. Useful for naming files."""
    #TODO: should be able to handle a time range
    time = parse_time(t)
    result = time.strftime("%Y%m%d_%H%M%S")
    return result
