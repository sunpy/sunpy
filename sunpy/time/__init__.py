"""Time related functionality"""
from datetime import datetime
from datetime import timedelta
from sunpy.time.timerange import TimeRange
from sunpy.time.julian import *

__all__ = ["parse_time", "day_of_year", "break_time"]
__all__ += julian.__all__
__all__ += timerange.__all__

def parse_time(time_string=None):
    """Given a time string will parse and return a datetime object.
    Similar to the anytim function in IDL.
    
    .. note :: If no input is given then returns the datetime for the 
    current time.
    
    Parameters
    ----------
    time_string : string
        
    Returns
    -------
    value : datetime

    See Also
    --------

    Examples
    --------
    >>> import sunpy.instr.goes as goes
    >>> goes.get_file(('2011/04/04', '2011/04/05'))
    
    .. todo:: add ability to parse tai (International Atomic Time seconds since 
    Jan 1, 1958)

    """
    if time_string is None:
        return datetime.now()
    elif isinstance(time_string, datetime):
        return time_string
    elif isinstance(time_string, tuple):
        return datetime(*time_string)
    elif isinstance(time_string, int) or isinstance(time_string, float):
        return datetime(1979, 1, 1) + timedelta(0, time_string)
    else:
        time_format_list = \
            ["%Y-%m-%dT%H:%M:%S.%f",    # Example 2007-05-04T21:08:12.999999
             "%Y/%m/%dT%H:%M:%S.%f",    # Example 2007/05/04T21:08:12.999999
             "%Y-%m-%dT%H:%M:%S.%fZ",   # Example 2007-05-04T21:08:12.999Z
             "%Y-%m-%dT%H:%M:%S",       # Example 2007-05-04T21:08:12
             "%Y%m%dT%H%M%S.%f",        # Example 20070504T210812.999999
             "%Y%m%dT%H%M%S",           # Example 20070504T210812
             "%Y/%m/%d %H:%M:%S",       # Example 2007/05/04 21:08:12
             "%Y/%m/%d %H:%M",          # Example 2007/05/04 21:08
             "%Y/%m/%d %H:%M:%S.%f",    # Example 2007/05/04 21:08:12.999999
             "%Y-%m-%d %H:%M:%S.%f",    # Example 2007-05-04 21:08:12.999999
             "%Y-%m-%d %H:%M:%S",       # Example 2007-05-04 21:08:12
             "%Y-%m-%d %H:%M",          # Example 2007-05-04 21:08
             "%Y-%b-%d %H:%M:%S",       # Example 2007-May-04 21:08:12
             "%Y-%b-%d %H:%M",          # Example 2007-May-04 21:08
             "%Y-%b-%d",                # Example 2007-May-04
             "%Y-%m-%d",                # Example 2007-05-04
             "%Y/%m/%d",                # Example 2007/05/04
             "%Y%m%d_%H%M%S"]           # Example 20070504_210812
        for time_format in time_format_list: 
            try: 
                return datetime.strptime(time_string, time_format)
            except ValueError:
                pass
    
        raise ValueError("%s is not a valid time string!" % time_string)
    
def day_of_year(t=None):
    """Returns the day of year."""
    SECONDS_IN_DAY = 60 * 60 * 24.0
    time = parse_time(t)
    time_diff = parse_time(t) - datetime(time.year, 1, 1, 0, 0, 0)
    result = time_diff.days + time_diff.seconds / SECONDS_IN_DAY
    return result

def break_time(t=None):
    """Given a time returns a string. Useful for naming files."""
    #TODO: should be able to handle a time range
    return parse_time(t).strftime("%Y%m%d_%H%M%S")
