"""
This module provies a collection of time handing functions.
"""
import re
import textwrap
from datetime import date, datetime
from functools import singledispatch

import numpy as np

import astropy.time
import astropy.units as u
from astropy.time import Time

# This is not called but imported to register it
from sunpy.time.utime import TimeUTime  # NOQA
from sunpy.util.decorators import add_common_docstring

__all__ = [
    'find_time', 'parse_time', 'is_time',
    'is_time_in_given_format', 'is_time_equal',
    'julian_centuries'
]

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
    "%Y-%m-%dT%H:%M:%S.%f",  # Example 2007-05-04T21:08:12.999999
    "%Y/%m/%dT%H:%M:%S.%f",  # Example 2007/05/04T21:08:12.999999
    "%Y-%m-%dT%H:%M:%S.%fZ",  # Example 2007-05-04T21:08:12.999Z
    "%Y-%m-%dT%H:%M:%S",  # Example 2007-05-04T21:08:12
    "%Y/%m/%dT%H:%M:%S",  # Example 2007/05/04T21:08:12
    "%Y%m%dT%H%M%S.%f",  # Example 20070504T210812.999999
    "%Y%m%dT%H%M%S",  # Example 20070504T210812
    "%Y/%m/%d %H:%M:%S",  # Example 2007/05/04 21:08:12
    "%Y/%m/%d %H:%M",  # Example 2007/05/04 21:08
    "%Y/%m/%d %H:%M:%S.%f",  # Example 2007/05/04 21:08:12.999999
    "%Y-%m-%d %H:%M:%S.%f",  # Example 2007-05-04 21:08:12.999999
    "%Y-%m-%d %H:%M:%S",  # Example 2007-05-04 21:08:12
    "%Y-%m-%d %H:%M",  # Example 2007-05-04 21:08
    "%Y-%b-%d %H:%M:%S",  # Example 2007-May-04 21:08:12
    "%Y-%b-%d %H:%M",  # Example 2007-May-04 21:08
    "%Y-%b-%d",  # Example 2007-May-04
    "%Y-%m-%d",  # Example 2007-05-04
    "%Y/%m/%d",  # Example 2007/05/04
    "%d-%b-%Y",  # Example 04-May-2007
    "%d-%b-%Y %H:%M:%S.%f",  # Example 04-May-2007 21:08:12.999999
    "%Y%m%d_%H%M%S",  # Example 20070504_210812
    "%Y:%j:%H:%M:%S",  # Example 2012:124:21:08:12
    "%Y:%j:%H:%M:%S.%f",  # Example 2012:124:21:08:12.999999
    "%Y%m%d%H%M%S",  # Example 20140101000001 (JSOC / VSO)
    "%Y.%m.%d_%H:%M:%S_TAI",  # Example 2016.05.04_21:08:12_TAI
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
        return inp, astropy.time.TimeDelta(0 * u.day)
    if hour == "24":
        if not all(
                _n_or_eq(_group_or_none(match, g, int), 00)
                for g in ["minute", "second", "microsecond"]):
            raise ValueError
        from_, to = match.span("hour")
        return inp[:from_] + "00" + inp[to:], astropy.time.TimeDelta(1 * u.day)
    return inp, astropy.time.TimeDelta(0 * u.day)


def find_time(string, format):
    """
    Return iterator of occurrences of date formatted with format in string.

    Currently supported format codes:
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


@singledispatch
def convert_time(time_string, format=None, **kwargs):
    # default case when no type matches
    return Time(time_string, format=format, **kwargs)


# Only register pandas if we can import pandas
try:
    import pandas

    @convert_time.register(pandas.Timestamp)
    def convert_time_pandasTimestamp(time_string, **kwargs):
        return Time(time_string.asm8)

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
    time_string = (time_string + (0, ) * 7)[:7]
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


@convert_time.register(list)
def convert_time_list(time_list, format=None, **kwargs):
    item = time_list[0]
    # If we have a list of strings, need to get the correct format from our
    # list of custom formats.
    if isinstance(item, str) and format is None:
        string_format = _get_time_fmt(item)
        return Time.strptime(time_list, string_format, **kwargs)

    # Otherwise return the default method
    return convert_time.dispatch(object)(time_list, format, **kwargs)


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
                ts, time_delta = _regex_parse_time(time_string, time_format)
            except TypeError:
                break
            if ts is None:
                continue
            return Time.strptime(ts, time_format, **kwargs) + time_delta
        except ValueError:
            pass

    # when no format matches, call default fucntion
    return convert_time.dispatch(object)(time_string, **kwargs)


def _get_time_fmt(time_string):
    """
    Try all the formats in TIME_FORMAT_LIST to work out which one applies to
    the time string.
    """
    for time_format in TIME_FORMAT_LIST:
        ts, _ = _regex_parse_time(time_string, time_format)
        if ts is not None:
            return time_format


def _variables_for_parse_time_docstring():
    ret = {}

    example_time = datetime(2017, 1, 1, 11, 10, 9)
    example_parse_time = [example_time.strftime(f) for f in TIME_FORMAT_LIST]
    example_parse_time = "\n      ".join(example_parse_time)
    ret['parse_time_formats'] = example_parse_time

    types = list(convert_time.registry.keys())
    types.remove(object)
    # Do Builtins
    types2 = [t.__qualname__ for t in types if t.__module__ == "builtins"]
    # # Do all the non-special ones where we take the package name and the class
    types2 += [t.__module__.split(".")[0] + "." +
               t.__qualname__ for t in types if not t.__module__.startswith(("builtins", "astropy"))]
    # Special case astropy.time where we need the subpackage
    types2 += ["astropy.time." +
               t.__qualname__ for t in types if t.__module__.startswith("astropy.time")]
    parse_time_types = str(types2)[1:-1].replace("'", "`")
    ret['parse_time_types'] = parse_time_types
    ret['parse_time_desc'] = """
                             Any time input, will be passed into `~sunpy.time.parse_time`.
                             """
    ret['astropy_time_formats'] = textwrap.fill(str(list(astropy.time.Time.FORMATS.keys())),
                                                subsequent_indent=' '*10)

    return ret


@add_common_docstring(**_variables_for_parse_time_docstring())
def parse_time(time_string, *, format=None, **kwargs):
    """
    Takes a time input and will parse and return a `astropy.time.Time`.

    Parameters
    ----------
    time_string : {parse_time_types}
        Time to parse.
    format : `str`, optional
        Specifies the format user has provided the time_string in.
        We support the same formats of `astropy.time.Time`.

        The allowed values for ``format`` can be listed with::

          >>> list(astropy.time.Time.FORMATS)
          {astropy_time_formats}

    Returns
    -------
    `astropy.time.Time`
        `~astropy.time.Time` corresponding to input time string.

    Examples
    --------
    >>> import sunpy.time
    >>> sunpy.time.parse_time('2012/08/01')
    <Time object: scale='utc' format='isot' value=2012-08-01T00:00:00.000>
    >>> sunpy.time.parse_time('2016.05.04_21:08:12_TAI')
    <Time object: scale='tai' format='isot' value=2016-05-04T21:08:12.000>

    Notes
    -----
    Additional keyword arguments can be passed to `astropy.time.Time`

    The list of time formats are show by the following examples::

      {parse_time_formats}
    """
    if isinstance(time_string, str) and time_string == 'now':
        rt = Time.now()
    else:
        rt = convert_time(time_string, format=format, **kwargs)

    return rt


def is_time(time_string, time_format=None):
    """
    Returns true if the input is a valid date/time representation.

    Parameters
    ----------
    time_string : `int`, `float`, ``time_string``, `datetime`
        Date to parse.
    time_format : basestring, utime, `datetime`, optional
        Specifies the format user has provided the ``time_string`` in.

    Returns
    -------
    `bool`
        `True` if can be parsed by `sunpy.time.parse_time`.

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


def is_time_in_given_format(time_string, time_format):
    """
    Tests whether a time string is formatted according to the given time
    format.

    Parameters
    ----------
    time_string : `str`
        Date to parse.
    time_format : basestring, utime, `datetime`, optional
        Specifies the format user has provided the ``time_string`` in.

    Returns
    -------
    `bool`
        `True` if can it if follows the ``time_format``.
    """
    try:
        datetime.strptime(time_string, time_format)
        return True
    except ValueError:
        return False


def julian_centuries(t='now'):
    """
    Returns the number of Julian centuries since J1900.0 (noon on 1900 January
    0).
    """
    DAYS_IN_JULIAN_CENTURY = 36525.0

    # J1900.0 is 2415021.0
    return (parse_time(t).jd - 2415020.0) / DAYS_IN_JULIAN_CENTURY
