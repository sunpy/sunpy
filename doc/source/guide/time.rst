.. _time-in-sunpy:

=============
Time in SunPy
=============

Working with times and time ranges is a standard task in solar data analysis as such
SunPy strives to provide convenient and easy methods to do the simple stuff. Python
already provides an object for a time or date through `datetime.datetime`.
SunPy builds upon its functionality.

.. _parse-time:

1. Parsing Times
----------------

Solar data is associated with a number of different time formats. SunPy provides a simple
parsing function which can deal with most every format that a user may encounter. Called
`sunpy.time.parse_time()`, this function takes a string as input and returns a datetime object.
Here are few examples of formats which `sunpy.time.parse_time()` accepts: ::

    >>> from sunpy.time import parse_time
    >>> parse_time('2007-05-04T21:08:12')   # doctest: +SKIP
    >>> parse_time('2007/05/04T21:08:12')   # doctest: +SKIP
    >>> parse_time('20070504T210812')   # doctest: +SKIP
    >>> parse_time('2007-May-04 21:08:12')   # doctest: +SKIP
    >>> parse_time('20070504_210812')   # doctest: +SKIP
    datetime.datetime(2007, 5, 4, 21, 8, 12)

Each of the above returns the same datetime object ``datetime.datetime(2007,
5, 4, 21, 8, 12)``. One of the most standard time formats used in solar
physics is the number of seconds since 1979 January 01. The parse_time
function also accepts this as input, e.g.: ::

    >>> parse_time(894316092.00000000)
    datetime.datetime(2007, 5, 4, 21, 8, 12)


All SunPy functions which require
time as an input sanitize the input using parse_time.

2. Time Ranges
--------------

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
The duration must be provided as a `datetime.timedelta` object or
time-equivalent `astropy.units.Quantity`
example: ::

    >>> import astropy.units as u
    >>> time_range = TimeRange('2010/03/04 00:10', 400 * u.second)

or: ::

    >>> from datetime import timedelta
    >>> time_range = TimeRange('2010/03/04 00:10', timedelta(0, 400))

The time range objects provides a number of useful functions. For example, you can easily
get the time at the center of your interval or the length of your interval in minutes
or days or seconds: ::

    >>> time_range.center
    datetime.datetime(2010, 3, 4, 0, 13, 20)
    >>> time_range.minutes
    <Quantity 6.666666666666667 min>
    >>> time_range.days
    <Quantity 0.004629629629629629 d>
    >>> time_range.seconds
    <Quantity 400.0 s>

It also makes it easy to create new time ranges. The functions next() and previous()
do an inplace update to the object by either adding or subtracting the same time interval
. This could be useful if you need to step through a number of time ranges. For example,
if you needed time ranges that spanned 30 minutes over a period of 4 hours you could do: ::

    >>> for a in range(8):
    ...     print(time_range.next())
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
