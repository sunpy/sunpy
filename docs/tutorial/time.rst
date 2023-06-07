.. _sunpy-tutorial-times:

*****
Times
*****

In this section of the tutorial, you will learn how times are represented and manipulated in sunpy.
sunpy makes extensive use of the `astropy.time` module for this task.
By the end of this section of the tutorial, you will learn how to parse times from different formats as well as construct and inspect time ranges.

.. _sunpy-tutorial-times-parse-time:

Parsing Times
=============

Solar data is associated with a number of different time formats.
To handle all these formats, sunpy has :meth:`sunpy.time.parse_time` that accepts a variety of inputs, and returns a consistent `~astropy.time.Time` object.
You might have come across another way of storing time that's built into Python itself, `datetime.datetime`.
`~datetime.datetime` does not provide support for common time formats used in solar physics or leap seconds, hence the use of `astropy.time.Time` throughout sunpy.

Here's a few examples of using `~sunpy.time.parse_time` to create `~astropy.time.Time` objects:

.. code-block:: python

    >>> from sunpy.time import parse_time

    >>> parse_time('2007-05-04T21:08:12')
    <Time object: scale='utc' format='isot' value=2007-05-04T21:08:12.000>
    >>> parse_time(894316092.00000000, format='utime')
    <Time object: scale='utc' format='utime' value=894316092.0>

See the documentation for `sunpy.time.parse_time` for a full list of allowed arguments.

Time Ranges
===========

Another standard task in data analysis is dealing with pairs of times or time ranges.
To deal with time ranges sunpy provides the `sunpy.time.TimeRange` object.
A `sunpy.time.TimeRange` object can be created by providing a start time and an end time:

.. code-block:: python

    >>> from sunpy.time import TimeRange

    >>> time_range = TimeRange('2010/03/04 00:10', '2010/03/04 00:20')

`~sunpy.time.TimeRange` makes use of :meth:`sunpy.time.parse_time` so it can accept a wide variety of time formats.
Alternatively, you can specify a start time and a duration:

.. code-block:: python

    >>> import astropy.units as u

    >>> time_range = TimeRange('2010/03/04 00:10', 400 * u.second)

The time range objects provides a number of useful functions.
For example, you can easily get the time at the center of your interval or the length of the interval:

.. code-block:: python

    >>> time_range.center
    <Time object: scale='utc' format='isot' value=2010-03-04T00:13:20.000>
    >>> time_range.seconds
    <Quantity 400. s>

A time range can also be easily split into sub-intervals of equal length, for example to
split a TimeRange object into two new TimeRange objects:

.. code-block:: python

    >>> time_range.split(2)
    [   <sunpy.time.timerange.TimeRange object ...>
        Start: 2010-03-04 00:10:00
        End:   2010-03-04 00:13:20
        Center:2010-03-04 00:11:40
        Duration:0.002314814814814825 days or
               0.0555555555555558 hours or
               3.333333333333348 minutes or
               200.00000000000088 seconds
    ,    <sunpy.time.timerange.TimeRange object ...>
        Start: 2010-03-04 00:13:20
        End:   2010-03-04 00:16:40
        Center:2010-03-04 00:15:00
        Duration:0.002314814814814825 days or
               0.0555555555555558 hours or
               3.333333333333348 minutes or
               200.00000000000088 seconds
    ]

Check out the code reference for the `sunpy.time.TimeRange` object for more information.
