import re
from datetime import datetime, date
from functools import singledispatch

import numpy as np

import astropy.time
from sunpy.time import Time
import astropy.units as u

from sunpy.time.utime import TimeUTime  # noqa: F401

__all__ = ['find_time', 'parse_time', 'is_time',
           'day_of_year', 'break_time', 'get_day', 'is_time_in_given_format',
           'is_time_equal']

# Mapping of time format codes to regular expressions.
REGEX = {
    '%Y': r'(?P<year>\d{4})',
    '%j': r'(?P<dayofyear>\d{3})',
    '%m': r'(?P<month>\d{1,2})',
    '%d': r'(?P<day>\d{1,2})',
    '%H': r'(?P<hour>\d{1,2})',
    '%M': r'(?P<minute>\d{1,2})',
    '%S': r'(?P<second>\d{1,2})',
    '%f': r'(?P<microsecond>\d+)',
    '%b': r'(?P<month_str>[a-zA-Z]+)',
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
    "%d-%b-%Y %H:%M:%S.%f",    # Example 04-May-2007 21:08:12.999999
    "%Y%m%d_%H%M%S",           # Example 20070504_210812
    "%Y:%j:%H:%M:%S",          # Example 2012:124:21:08:12
    "%Y:%j:%H:%M:%S.%f",       # Example 2012:124:21:08:12.999999
    "%Y%m%d%H%M%S",            # Example 20140101000001 (JSOC / VSO)
    "%Y.%m.%d_%H:%M:%S_TAI",   # Example 2016.05.04_21:08:12_TAI
]


def is_time_equal(t1, t2):
    """
    Work around for https://github.com/astropy/astropy/issues/6970.
    Remove the usage of this function once the fix is in place.
    """
    if abs(t1 - t2) < 1 * u.nanosecond:
        return True
    return False


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
    for key, value in REGEX.items():
        format = format.replace(key, value)
    match = re.match(format, inp)
    if match is None:
        return None, None
    try:
        hour = match.group("hour")
    except IndexError:
        return inp, astropy.time.TimeDelta(0*u.day)
    if hour == "24":
        if not all(
                   _n_or_eq(_group_or_none(match, g, int), 00)
                   for g in ["minute", "second", "microsecond"]
                  ):
            raise ValueError
        from_, to = match.span("hour")
        return inp[:from_] + "00" + inp[to:], astropy.time.TimeDelta(1*u.day)
    return inp, astropy.time.TimeDelta(0*u.day)


def find_time(string, format):
    """
    Return iterator of occurrences of date formatted with format
    in string.

    Currently supported format codes: TODO: ADD THIS
    """
    re_format = format
    for key, value in REGEX.items():
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


find_time.__doc__ += ', '.join(list(REGEX.keys()))


def _iter_empty(iter):
    try:
        next(iter)
    except StopIteration:
        return True
    return False


def _astropy_time(time):
    """
    Return an `~astropy.time.Time` instance, running it through `~sunpy.time.parse_time` if needed
    """
    return time if isinstance(time, astropy.time.Time) else astropy.time.Time(parse_time(time))


@singledispatch
def convert_time(time_string, format=None, **kwargs):
    # default case when no type matches
    return Time(time_string, format=format, **kwargs)


# Only register pandas if we can import pandas
try:
    import pandas

    @convert_time.register(pandas.Timestamp)
    def convert_time_pandasTimestamp(time_string, **kwargs):
        return Time(time_string.to_pydatetime())

    @convert_time.register(pandas.Series)
    def convert_time_pandasSeries(time_string, **kwargs):
        return Time(time_string.tolist(), **kwargs)

    @convert_time.register(pandas.DatetimeIndex)
    def convert_time_pandasDatetimeIndex(time_string, **kwargs):
        return Time(time_string.tolist(), **kwargs)

except ImportError:
    pass


@convert_time.register(datetime)
def convert_time_datetime(time_string, **kwargs):
    return Time(time_string, **kwargs)


@convert_time.register(date)
def convert_time_date(time_string, **kwargs):
    return Time(time_string.isoformat(), **kwargs)


@convert_time.register(tuple)
def convert_time_tuple(time_string, **kwargs):
    # Make sure there are enough values to unpack
    time_string = (time_string + (0,)*7)[:7]
    return Time('{}-{}-{}T{}:{}:{}.{:06}'.format(*time_string), **kwargs)


@convert_time.register(np.datetime64)
def convert_time_npdatetime64(time_string, **kwargs):
    return Time(str(time_string.astype('M8[ns]')), **kwargs)


@convert_time.register(np.ndarray)
def convert_time_npndarray(time_string, **kwargs):
    if 'datetime64' in str(time_string.dtype):
        return Time([str(dt.astype('M8[ns]')) for dt in time_string], **kwargs)
    else:
        return convert_time.dispatch(object)(time_string, **kwargs)


@convert_time.register(astropy.time.Time)
def convert_time_astropy(time_string, **kwargs):
    return time_string


@convert_time.register(str)
def convert_time_str(time_string, **kwargs):
    # remove trailing zeros and the final dot to allow any
    # number of zeros. This solves issue #289
    if '.' in time_string:
        time_string = time_string.rstrip("0").rstrip(".")

    if 'TAI' in time_string:
        kwargs['scale'] = 'tai'

    for time_format in TIME_FORMAT_LIST:
        try:
            try:
                ts, time_delta = _regex_parse_time(time_string,
                                                   time_format)
            except TypeError:
                break
            if ts is None:
                continue
            return Time.strptime(ts, time_format, **kwargs) + time_delta
        except ValueError:
            pass

    # when no format matches, call default fucntion
    return convert_time.dispatch(object)(time_string, **kwargs)


def parse_time(time_string, *, format=None, **kwargs):
    """Given a time string will parse and return a `astropy.time.Time` object.
    Similar to the anytim function in IDL.
    utime -- Time since epoch 1 Jan 1979

    Parameters
    ----------
    time_string : [ int, float, string, datetime, astropy.time.Time,
                    numpy.datetime64, pandas.Timestamp ]
        Time to parse.

    format : [ 'jd', 'mjd', 'decimalyear', 'unix', 'cxcsec', 'gps',
               'plot_date', 'datetime', 'iso', 'isot', 'yday', 'fits',
               'byear', 'jyear', 'byear_str', 'jyear_str', 'utime']

        Specifies the format user has provided the time_string in.
        Same as format of `astropy.time.Time`

    kwargs : dict
        Additional keyword arguments that can be passed to `astropy.time.Time`

    Returns
    -------
    out : Time
        `astropy.time.Time` corresponding to input time string

    Examples
    --------
    >>> import sunpy.time
    >>> sunpy.time.parse_time('2012/08/01')
    <Time object: scale='utc' format='isot' value=2012-08-01T00:00:00.000>
    >>> sunpy.time.parse_time('2016.05.04_21:08:12_TAI')
    <Time object: scale='tai' format='isot' value=2016-05-04T21:08:12.000>
    """
    if time_string is 'now':
        rt = Time.now()
    else:
        rt = convert_time(time_string, format=format, **kwargs)

    return rt


def is_time(time_string, time_format=None):
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

    Examples
    --------
    >>> import sunpy.time
    >>> sunpy.time.is_time('2012/08/01')
    True
    """
    if time_string is None:
        return False
    elif isinstance(time_string, Time):
        return True

    try:
        parse_time(time_string, format=time_format)
    except ValueError:
        return False
    else:
        return True


def day_of_year(time_string):
    """
    Returns the (fractional) day of year.

    Note: This function takes into account leap seconds.

    Parameters
    ----------
    time_string : str
        A parse_time compatible string

    Returns
    -------
    out : float
        The fractional day of year (where Jan 1st is 1).

    Examples
    --------
    >>> import sunpy.time
    >>> sunpy.time.day_of_year('2012/01/01')
    1.0
    >>> sunpy.time.day_of_year('2012/08/01')
    214.00001157407408
    >>> sunpy.time.day_of_year('2005-08-04T00:18:02.000Z')
    216.01252314814815

    """
    SECONDS_IN_DAY = 60 * 60 * 24.0
    time = parse_time(time_string)
    time_diff = time - Time.strptime(time.strftime('%Y'), '%Y')
    return time_diff.jd + 1


def break_time(t='now', time_format=None):
    """Given a time returns a string. Useful for naming files."""
    # TODO: should be able to handle a time range
    return parse_time(t, format=time_format).strftime("%Y%m%d_%H%M%S")


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
