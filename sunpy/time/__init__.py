"""Time related functionality"""
from datetime import datetime
from datetime import timedelta
from . util import *

__all__ = []
__all__ += util.__all__

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
    
    Reference
    ---------
    | 
    
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
