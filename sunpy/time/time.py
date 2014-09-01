import re
from datetime import datetime
from datetime import timedelta

__all__ = ['find_time', 'extract_time', 'parse_time', 'is_time', 'day_of_year', 'break_time', 'get_day', 'is_time_in_given_format']

# Mapping of time format codes to regular expressions.
REGEX = {
    '%Y': '(?P<year>\d{4})',
    '%j': '(?P<dayofyear>\d{3})',
    '%m': '(?P<month>\d{1,2})',
    '%d': '(?P<day>\d{1,2})',
    '%H': '(?P<hour>\d{1,2})',
    '%M': '(?P<minute>\d{1,2})',
    '%S': '(?P<second>\d{1,2})',
    '%f': '(?P<microsecond>\d+)',
    '%b': '(?P<month_str>[a-zA-Z]+)',
}

TIME_FORMAT_LIST = [
    "%Y-%m-%dT%H:%M:%S.%f",    # Example 2007-05-04T21:08:12.999999
    "%Y/%m/%dT%H:%M:%S.%f",    # Example 2007/05/04T21:08:12.999999
    "%Y-%m-%dT%H:%M:%S.%fZ",   # Example 2007-05-04T21:08:12.999Z
    "%Y-%m-%dT%H:%M:%S",       # Example 2007-05-04T21:08:12
    "%Y/%m/%dT%H:%M:%S",       # Example 2007/05/04T21:08:12
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
    "%d-%b-%Y",                # Example 04-May-2007
    "%Y%m%d_%H%M%S",           # Example 20070504_210812
    "%Y:%j:%H:%M:%S",          # Example 2012:124:21:08:12
    "%Y:%j:%H:%M:%S.%f",       # Example 2012:124:21:08:12.999999
]


def _group_or_none(match, group, fun):
    try:
        ret = match.group(group)
    except IndexError:
        return None
    else:
        return fun(ret)


def _n_or_eq(a, b):
    return a is None or a == b


def _regex_parse_time(inp, format):
    # Parser for finding out the minute value so we can adjust the string
    # from 24:00:00 to 00:00:00 the next day because strptime does not
    # understand the former.
    for key, value in REGEX.iteritems():
        format = format.replace(key, value)
    match = re.match(format, inp)
    if match is None:
        return None, None
    try:
        hour = match.group("hour")
    except IndexError:
        return inp, timedelta(days=0)
    if match.group("hour") == "24":
        if not all(_n_or_eq(_group_or_none(match, g, int), 00)
            for g in ["minute", "second", "microsecond"]
        ):
            raise ValueError
        from_, to = match.span("hour")
        return inp[:from_] + "00" + inp[to:], timedelta(days=1)
    return inp, timedelta(days=0)


def find_time(string, format):
    """ Return iterator of occurences of date formatted with format
    in string. Currently supported format codes: """
    re_format = format
    for key, value in REGEX.iteritems():
        re_format = re_format.replace(key, value)
    matches = re.finditer(re_format, string)
    for match in matches:
        try:
            matchstr = string[slice(*match.span())]
            dt = datetime.strptime(matchstr, format)
        except ValueError:
            continue
        else:
            yield dt


find_time.__doc__ += ', '.join(REGEX.keys())


def _iter_empty(iter):
    try:
        iter.next()
    except StopIteration:
        return True
    return False


def extract_time(string):
    """ Find subset of string that corresponds to a datetime and return
    its value as a a datetime. If more than one or none is found, raise
    ValueError. """
    matched = None
    bestmatch = None
    for time_format in TIME_FORMAT_LIST:
        found = find_time(string, time_format)
        try:
            match = found.next()
        except StopIteration:
            continue
        else:
            if matched is not None:
                if time_format.startswith(matched):
                    # Already matched is a substring of the one just matched.
                    matched = time_format
                    bestmatch = match
                elif not matched.startswith(time_format):
                    # If just matched is substring of time_format, just ignore
                    # just matched.
                    raise ValueError("Ambiguous string")
            else:
                matched = time_format
                bestmatch = match
            if not _iter_empty(found):
                raise ValueError("Ambiguous string")
    if not matched:
        raise ValueError("Time not found")
    return bestmatch


def parse_time(time_string, time_format=''):
    """Given a time string will parse and return a datetime object.
    Similar to the anytim function in IDL.
    utime -- Time since epoch 1 Jan 1979

    Parameters
    ----------
    time_string : [ int, float, time_string, datetime ]
        Date to parse which can be either time_string, int, datetime object.
    time_format : [ basestring, utime, datetime ]
        Specifies the format user has provided the time_string in.

    Returns
    -------
    out : datetime
        DateTime corresponding to input date string

    Note:
    If time_string is an instance of float, then it is assumed to be in utime format.

    Examples
    --------
    >>> sunpy.time.parse_time('2012/08/01')
    >>> sunpy.time.parse_time('2005-08-04T00:01:02.000Z')

    Todo:
    Add ability to parse tai (International Atomic Time seconds since
    Jan 1, 1958)
    """
    if isinstance(time_string, datetime) or time_format == 'datetime':
        return time_string
    elif isinstance(time_string, tuple):
        return datetime(*time_string)
    elif time_format == 'utime' or  isinstance(time_string, (int, float))  :
        return datetime(1979, 1, 1) + timedelta(0, time_string)
    elif time_string is 'now':
        return datetime.utcnow()
    else:
        # remove trailing zeros and the final dot to allow any
        # number of zeros. This solves issue #289
        if '.' in time_string:
            time_string = time_string.rstrip("0").rstrip(".")
        for time_format in TIME_FORMAT_LIST:
            try:
                try:
                    ts, time_delta = _regex_parse_time(time_string,
                                                       time_format)
                except TypeError:
                    break
                if ts is None:
                    continue
                return datetime.strptime(ts, time_format) + time_delta
            except ValueError:
                pass
        raise ValueError("%s is not a valid time string!" % time_string)


def is_time(time_string, time_format=''):
    """
    Returns true if the input is a valid date/time representation

    Parameters
    ----------
    time_string : [ int, float, time_string, datetime ]
        Date to parse which can be either time_string, int, datetime object.
    time_format : [ basestring, utime, datetime ]
	Specifies the format user has provided the time_string in.
   
    Returns
    -------
    out : bool
        True if can be parsed by parse_time

    Note:
    If time_string is an instance of float, then it is assumed to be in unix time format.
    Examples
    --------
    >>> sunpy.time.parse_time('2012/08/01')
    >>> sunpy.time.parse_time('2005-08-04T00:01:02.000Z')

    .. todo::

        add ability to parse tai (International Atomic Time seconds
        since Jan 1, 1958)

    """
    if time_string is None:
        return False
    elif isinstance(time_string, datetime):
        return True

    try:
        parse_time(time_string,time_format)
    except ValueError:
        return False
    else:
        return True


def day_of_year(time_string):
    """Returns the (fractional) day of year.

    Parameters
    ----------
    time_string : string
        A parse_time compatible string

    Returns
    -------
    out : float
        The fractional day of year (where Jan 1st is 1).

    Examples
    --------
    >>> sunpy.time.day_of_year('2012/01/01')
    1.00
    >>> sunpy.time.day_of_year('2012/08/01')
    214.00
    >>> sunpy.time.day_of_year('2005-08-04T00:18:02.000Z')
    216.01252314814815

    """
    SECONDS_IN_DAY = 60 * 60 * 24.0
    time = parse_time(time_string)
    time_diff = time - datetime(time.year, 1, 1, 0, 0, 0)
    return time_diff.days + time_diff.seconds / SECONDS_IN_DAY + 1

def break_time(t='now', time_format=''):
    """Given a time returns a string. Useful for naming files."""
    #TODO: should be able to handle a time range
    return parse_time(t, time_format).strftime("%Y%m%d_%H%M%S")

def get_day(dt):
    """ Return datetime for the beginning of the day of given datetime. """
    return datetime(dt.year, dt.month, dt.day)

def is_time_in_given_format(time_string, time_format):
    """Tests whether a time string is formatted according to the given time
    format."""
    try:
        datetime.strptime(time_string, time_format)
        return True
    except ValueError:
        return False
