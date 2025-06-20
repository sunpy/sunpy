"""
This module provides a collection of time handing functions.
"""
import re
import textwrap
import contextlib
from datetime import date, datetime
from functools import singledispatch

import numpy as np

import astropy.time
import astropy.units as u
from astropy.time import Time, TimeDelta
from astropy.time import conf as astropy_time_conf
from astropy.utils.introspection import minversion

from sunpy import log
# This is not called but imported to register time formats
from sunpy.time.timeformats import *  # NOQA
from sunpy.util.decorators import add_common_docstring

__all__ = [
    'find_time', 'parse_time', 'is_time',
    'is_time_in_given_format', 'is_time_equal',
    'julian_centuries'
]

# Mapping of time format codes to regular expressions.
REGEX = {
    '.': r'\.',
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
# DO NOT SORT THIS LIST.
# The string parsing is dependent on the specific order within this list.
TIME_FORMAT_LIST = [
    "%Y-%m-%dT%H:%M:%S.%f",  # Example 2007-05-04T21:08:12.999999
    "%Y/%m/%dT%H:%M:%S.%f",  # Example 2007/05/04T21:08:12.999999
    "%Y-%m-%dT%H:%M:%S.%fZ",  # Example 2007-05-04T21:08:12.999Z
    "%Y-%m-%dT%H:%M:%S",  # Example 2007-05-04T21:08:12
    "%Y/%m/%dT%H:%M:%S",  # Example 2007/05/04T21:08:12
    "%Y-%m-%dT%H:%M:%SZ",  # Example 2007-05-04T21:08:12Z
    "%Y%m%dT%H%M%S.%f",  # Example 20070504T210812.999999
    "%Y%m%dT%H%M",  # Example 20070504T2108 , Should precede "%Y%m%dT%H%M%S".
    "%Y%m%dT%H%M%S",  # Example 20070504T210812
    "%Y/%m/%d %H:%M:%S",  # Example 2007/05/04 21:08:12
    "%Y/%m/%d %H:%M",  # Example 2007/05/04 21:08
    "%Y/%m/%d %H:%M:%S.%f",  # Example 2007/05/04 21:08:12.999999
    "%Y-%m-%d %H:%M:%S.%f",  # Example 2007-05-04 21:08:12.999999
    "%Y-%m-%d %H:%M:%S",  # Example 2007-05-04 21:08:12
    "%Y-%m-%d %H:%M",  # Example 2007-05-04 21:08
    "%Y-%b-%d %H:%M:%S.%f",  # Example 2007-May-04 21:08:12.999999
    "%Y-%b-%d %H:%M:%S",  # Example 2007-May-04 21:08:12
    "%Y-%b-%d %H:%M",  # Example 2007-May-04 21:08
    "%Y-%b-%d",  # Example 2007-May-04
    "%Y-%m-%d",  # Example 2007-05-04
    "%Y/%m/%d",  # Example 2007/05/04
    "%d-%b-%Y",  # Example 04-May-2007
    "%d-%b-%Y %H:%M:%S",  # Example 04-May-2007 21:08:12
    "%d-%b-%Y %H:%M:%S.%f",  # Example 04-May-2007 21:08:12.999999
    "%Y%m%d_%H%M",  # Example 20070504_2108 , Should precede "%Y%m%d_%H%M%S".
    "%Y%m%d_%H%M%S",  # Example 20070504_210812
    "%Y:%j:%H:%M:%S",  # Example 2012:124:21:08:12
    "%Y:%j:%H:%M:%S.%f",  # Example 2012:124:21:08:12.999999
    "%Y%m%d%H%M" ,   # Example 201401041205 , Should precede "%Y%m%d%H%M%S".
    "%Y%m%d%H%M%S",  # Example 20140101000001 (JSOC/VSO Export/Downloads)
    "%Y.%m.%d_%H:%M:%S_TAI",  # Example 2016.05.04_21:08:12_TAI - JSOC
    "%Y.%m.%d_%H:%M:%S.%f_TAI",  # Example 2019.09.15_00:00:02.898_TAI - JSOC
    "%Y.%m.%d_%H:%M:%S_UTC",  # Example 2016.05.04_21:08:12_UTC - JSOC
    "%Y.%m.%d_%H:%M:%S",  # Example 2016.05.04_21:08:12 - JSOC
    "%Y/%m/%dT%H:%M",  # Example 2007/05/04T21:08
]

_ONE_DAY_TIMEDELTA = TimeDelta(1 * u.day)


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
    match = re.match(f"{format}$", inp)
    if match is None:
        return None, None

    found_groups = match.groupdict()

    # Special handling to strip any excess zeros beyond the six digits of the microsecond field
    if "microsecond" in found_groups and re.match(r"\d{6}0+", match.group("microsecond")):
        from_, to = match.span("microsecond")
        inp = inp[:from_] + match.group("microsecond")[:6] + inp[to:]
        match = re.match(format, inp)

    # Special handling to add a day if the hour is 24 and the minute/second/microsecond are all zero
    add_one_day = False
    if "hour" in found_groups and match.group("hour") == "24":
        if not all(
                _n_or_eq(_group_or_none(match, g, int), 00)
                for g in ["minute", "second", "microsecond"]):
            raise ValueError
        from_, to = match.span("hour")
        inp = inp[:from_] + "00" + inp[to:]
        add_one_day = True
        match = re.match(format, inp)

    return inp, add_one_day


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


@convert_time.register(str)
@convert_time.register(list)
def convert_time_str(time_string, **kwargs):
    is_single_string = isinstance(time_string, str)
    is_string_list = isinstance(time_string, list) and isinstance(time_string[0], str)

    if is_single_string or is_string_list:
        first_item = time_string[0] if is_string_list else time_string
        if 'TAI' in first_item:
            kwargs['scale'] = 'tai'

        for time_format in TIME_FORMAT_LIST:
            try:
                try:
                    ts, add_one_day = _regex_parse_time(first_item, time_format)
                except TypeError:
                    break
                if ts is None:
                    continue
                if is_single_string:
                    t = Time.strptime(ts, time_format, **kwargs)
                    if add_one_day:
                        t += _ONE_DAY_TIMEDELTA
                else:
                    # For a list of strings, we do not try to correct 24:00:00
                    t = Time.strptime(time_string, time_format, **kwargs)
                return t
            except ValueError:
                pass

        log.debug("No matching sunpy format found for %s, so falling back to astropy formats", first_item)

    # If the string format does not match one of ours, we need to protect against a bad interaction
    # between astropy's C fast parser and numpy>=2.3
    # https://github.com/astropy/astropy/issues/18254
    # TODO: once the bug is fixed, skip this protection for fixed versions of astropy/numpy
    if is_string_list and len(time_string) > 500 and astropy_time_conf.use_fast_parser == "True" and minversion(np, "2.3.0"):
        log.debug("Disabling astropy's fast time parsing due to a known bug with long lists")
        ctx = astropy_time_conf.set_temp("use_fast_parser", "False")
    else:
        ctx = contextlib.nullcontext()

    with ctx:
        # when no format matches, call default function
        return convert_time.dispatch(object)(time_string, **kwargs)


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
    try:
        # Need to try importing cdflib, as if it is present it will register
        # extra formats with time
        import cdflib  # NOQA
    except Exception:
        pass
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
        We support the same formats as `astropy.time.Time`, which are::

          >>> list(astropy.time.Time.FORMATS)
          {astropy_time_formats}

    **kwargs :
        Additional keyword arguments are passed to `astropy.time.Time`

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
