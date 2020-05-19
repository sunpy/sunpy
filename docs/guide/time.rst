.. _time-in-sunpy:

*************
Time in SunPy
*************

Working with times and time ranges is a standard task in solar data analysis and as such
SunPy strives to provide convenient and easy methods to do the simple stuff. Python
already provides an object for a time or date through `datetime.datetime`.
However, `datetime.datetime` does not provide support for common time formats used in
solar physics nor leap seconds. To alleviate this, we use `astropy.time.Time` internally
which allows us to provide a superior user experience.

.. _parse-time:

1. Parsing Times
================

Solar data is associated with a number of different time formats. SunPy provides a simple
parsing function which can deal with most every format that a user may encounter. Called
`sunpy.time.parse_time()`, this function can take a variety of inputs.

Strings
-------

The most commonly used are strings and we support a selection of formats
which are matched using regrex. We currently support the following style of string formats::

    "2007-05-04T21:08:12.999999"
    "2007/05/04T21:08:12.999999"
    "2007-05-04T21:08:12.999Z"
    "2007-05-04T21:08:12"
    "2007/05/04T21:08:12"
    "20070504T210812.999999"
    "20070504T210812"
    "2007/05/04 21:08:12"
    "2007/05/04 21:08"
    "2007/05/04 21:08:12.999999"
    "2007-05-04 21:08:12.999999"
    "2007-05-04 21:08:12"
    "2007-05-04 21:08"
    "2007-May-04 21:08:12"
    "2007-May-04 21:08"
    "2007-May-04"
    "2007-05-04"
    "2007/05/04"
    "04-May-2007"
    "04-May-2007 21:08:12.999999"
    "20070504_210812"
    "2012:124:21:08:12"
    "2012:124:21:08:12.999999"
    "20140101000001"
    "2016.05.04_21:08:12_TAI"

If we pass some of these strings into `sunpy.time.parse_time()`::

    >>> from sunpy.time import parse_time
    >>> parse_time('2007-05-04T21:08:12')
    <Time object: scale='utc' format='isot' value=2007-05-04T21:08:12.000>
    >>> parse_time('2007/05/04T21:08:12')
    <Time object: scale='utc' format='isot' value=2007-05-04T21:08:12.000>
    >>> parse_time('20070504T210812')
    <Time object: scale='utc' format='isot' value=2007-05-04T21:08:12.000>
    >>> parse_time('2007-May-04 21:08:12')
    <Time object: scale='utc' format='isot' value=2007-05-04T21:08:12.000>
    >>> parse_time('20070504_210812')
    <Time object: scale='utc' format='isot' value=2007-05-04T21:08:12.000>

Each of the above returns the same `~astropy.time.Time` object ``<Time object: scale='utc' format='isot' value=2007-05-04T21:08:12.000>``.

We also support ``utime``, which is the amount of seconds from `1979-01-01 00:00:00 UTC`.
Same as Unix time but this starts 9 years later. The parse_time function also accepts this as input, e.g.::

    >>> parse_time(894316092.00000000, format='utime')
    <Time object: scale='utc' format='utime' value=894316092.0>

Other formats
-------------

`sunpy.time.parse_time()` understands more than just a string.
For example::

    >>> import datetime
    >>> parse_time(datetime.datetime(2007, 5, 4, 21, 8, 12))  # datetime.datetime
    <Time object: scale='utc' format='datetime' value=2007-05-04 21:08:12>
    >>> parse_time(datetime.date(2007, 5, 4))  # datetime.date
    <Time object: scale='utc' format='iso' value=2007-05-04 00:00:00.000>

    >>> parse_time((2007, 5, 4, 21, 8, 12))  # tuple of numbers
    <Time object: scale='utc' format='isot' value=2007-05-04T21:08:12.000>

    >>> import pandas
    >>> parse_time(pandas.Timestamp('2007-05-04T21:08:12'))  # pandas.Timestamp
    <Time object: scale='utc' format='datetime' value=2007-05-04 21:08:12>

    >>> time_ranges = [datetime.datetime(2007, 5, i) for i in range(1, 3)]
    >>> parse_time(pandas.Series(time_ranges))  # pandas.Series
    <Time object: scale='utc' format='datetime' value=[datetime.datetime(2007, 5, 1, 0, 0) datetime.datetime(2007, 5, 2, 0, 0)]>

    >>> parse_time(pandas.DatetimeIndex(time_ranges))  # pandas.DatetimeIndex
    <Time object: scale='utc' format='datetime' value=[datetime.datetime(2007, 5, 1, 0, 0) datetime.datetime(2007, 5, 2, 0, 0)]>

    >>> import numpy as np
    >>> parse_time(np.datetime64('2007-05-04T00'))  # np.datetime64
    <Time object: scale='utc' format='isot' value=2007-05-04T00:00:00.000>
    >>> parse_time(np.arange('2007-05-03', '2007-05-04', dtype='datetime64[D]'))  # np.ndarray
    <Time object: scale='utc' format='isot' value=['2007-05-03T00:00:00.000']>

`astropy.time.Time` API comparision
-----------------------------------

`sunpy.time.parse_time` is a wrapper around `astropy.time.Time`. The API is
nearly identical as `~astropy.time.Time` but supports more time input formats.
You can specify the format, scale, precision, location and other arguments just
as you would do with `~astropy.time.Time`. An example::

    >>> times = ['1999-01-01T00:00:00.123456789', '2010-01-01T00:00:00']
    >>> parse_time(times, format='isot', scale='tai')
    <Time object: scale='tai' format='isot' value=['1999-01-01T00:00:00.123' '2010-01-01T00:00:00.000']>

Please be aware that all SunPy functions which require time as an input sanitize the input using `~sunpy.time.parse_time`.

2. Time Ranges
==============

A very standard task in data analysis is to have to deal with pairs of times or time
ranges. This occurs very often with plotting or when searching for data. To deal with
time ranges SunPy provides the `sunpy.time.TimeRange` object. A TimeRange object can be created
very easily by providing it with two time strings, a start time and an end time: ::

    >>> from sunpy.time import TimeRange
    >>> time_range = TimeRange('2010/03/04 00:10', '2010/03/04 00:20')

You can also pass the start and end times as a tuple: ::

    >>> time_range = TimeRange(('2010/03/04 00:10', '2010/03/04 00:20'))

This object makes use of parse_time() so it can accept a wide variety of time formats.
A time range object can also be created by providing a start time and a duration.
The duration must be provided as a `~astropy.time.TimeDelta` or
time-equivalent `astropy.units.Quantity` or `datetime.timedelta` object
example: ::

    >>> import astropy.units as u
    >>> time_range = TimeRange('2010/03/04 00:10', 400 * u.second)

or: ::

    >>> import astropy.units as u
    >>> from astropy.time import TimeDelta
    >>> time_range = TimeRange('2010/03/04 00:10', TimeDelta(400 * u.second))

or: ::

    >>> from datetime import timedelta
    >>> time_range = TimeRange('2010/03/04 00:10', timedelta(0, 400))

The time range objects provides a number of useful functions. For example, you can easily
get the time at the center of your interval or the length of your interval in minutes
or days or seconds: ::

    >>> time_range.center
    <Time object: scale='utc' format='isot' value=2010-03-04T00:13:20.000>
    >>> time_range.minutes
    <Quantity 6.66666667 min>
    >>> time_range.days
    <Quantity 0.00462963 d>
    >>> time_range.seconds
    <Quantity 400. s>

It also makes it easy to create new time ranges. The functions next() and previous()
do an inplace update to the object by either adding or subtracting the same time interval
. This could be useful if you need to step through a number of time ranges. For example,
if you needed time ranges that spanned 30 minutes over a period of 4 hours you could do: ::

    >>> for a in range(8):
    ...     print(time_range.next())  # doctest: +IGNORE_OUTPUT
        Start: 2010-03-04 00:16:40
        End:   2010-03-04 00:23:20
        Center:2010-03-04 00:20:00
        Duration:0.00462962962963 days or
               0.111111111111 hours or
               6.66666666667 minutes or
               400.0 seconds
    <BLANKLINE>
        Start: 2010-03-04 00:23:20
        End:   2010-03-04 00:30:00
        Center:2010-03-04 00:26:40
        Duration:0.00462962962963 days or
               0.111111111111 hours or
               6.66666666667 minutes or
               400.0 seconds
    <BLANKLINE>
        Start: 2010-03-04 00:30:00
        End:   2010-03-04 00:36:40
        Center:2010-03-04 00:33:20
        Duration:0.00462962962963 days or
               0.111111111111 hours or
               6.66666666667 minutes or
               400.0 seconds
    <BLANKLINE>
        Start: 2010-03-04 00:36:40
        End:   2010-03-04 00:43:20
        Center:2010-03-04 00:40:00
        Duration:0.00462962962963 days or
               0.111111111111 hours or
               6.66666666667 minutes or
               400.0 seconds
    <BLANKLINE>
        Start: 2010-03-04 00:43:20
        End:   2010-03-04 00:50:00
        Center:2010-03-04 00:46:40
        Duration:0.00462962962963 days or
               0.111111111111 hours or
               6.66666666667 minutes or
               400.0 seconds
    <BLANKLINE>
        Start: 2010-03-04 00:50:00
        End:   2010-03-04 00:56:40
        Center:2010-03-04 00:53:20
        Duration:0.00462962962963 days or
               0.111111111111 hours or
               6.66666666667 minutes or
               400.0 seconds
    <BLANKLINE>
        Start: 2010-03-04 00:56:40
        End:   2010-03-04 01:03:20
        Center:2010-03-04 01:00:00
        Duration:0.00462962962963 days or
               0.111111111111 hours or
               6.66666666667 minutes or
               400.0 seconds
    <BLANKLINE>
        Start: 2010-03-04 01:03:20
        End:   2010-03-04 01:10:00
        Center:2010-03-04 01:06:40
        Duration:0.00462962962963 days or
               0.111111111111 hours or
               6.66666666667 minutes or
               400.0 seconds
    <BLANKLINE>

A time range can also be easily split into sub-intervals of equal length, for example to
split a TimeRange object into two new TimeRange objects: ::

    time_range.split(2)

Check out the code reference for the `sunpy.time.TimeRange` object for more information.
