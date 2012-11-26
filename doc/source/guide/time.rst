=============
Time in SunPy
=============

Working with times and time ranges is a standard task in solar data analysis as such
SunPy strives to provide convenient and easy methods to do the simple stuff. Python
already provides an object for a time or date through the `datetime object 
<http://docs.python.org/2/library/datetime.html>`_. SunPy builds upon its functionality.

1. Parsing Times
----------------

Solar data is associated with a number of different time formats. SunPy provides a simple
parsing function which can deal with most every format that a user may encounter. Called
parse_time(), this function takes a string as input and returns a datetime object.
Here are few examples of formats which parse_time() accepts: ::

    from sunpy.time import *
    parse_time('2007-05-04T21:08:12')
    parse_time('2007/05/04T21:08:12')
    parse_time('20070504T210812')
    parse_time('2007-May-04 21:08:12')
    parse_time('20070504_210812')

Each of the above returns the same datetime object. All SunPy functions which require 
time as an input sanitize the input using parse_time. You can also pass it a datetime
object directly and it will simply hand it right back to you. For users of IDL, 
this function is meant to the be the equivalent to SSW's anytim().

2. Time Ranges
--------------

A very standard task in data analysis is to have to deal with pairs of times or time 
ranges. This occurs very often with plotting or when searching for data. To deal with 
time ranges SunPy provides the TimeRange object. A time range object can be created
very easily by providing it with two time strings, a start time and an end time: ::

    time_range = TimeRange('2010/03/04 00:10', '2010/03/04 00:20')

You can also pass the start and end times as a tuple: ::

    time_range = TimeRange(('2010/03/04 00:10', '2010/03/04 00:20'))

This object makes use of parse_time() so it can accept a wide variety of time formats.
A time range object can also be created by providing a start time and a duration.
The duration must be provided as a `timedelta 
<http://docs.python.org/2/library/datetime.html#datetime.timedelta>`_ object or as a number of seconds so for
example: ::

    time_range = TimeRange('2010/03/04 00:10', 400)

or: ::

    from datetime import timedelta
    time_range = TimeRange('2010/03/04 00:10', timedelta(0,400))

The time range objects provides a number of useful functions. For example, you can easily
get the time at the center of your interval or the length of your interval in minutes 
or days or seconds: ::

    time_range.center()
    time_range.minutes()
    time_range.days()
    time_range.seconds()
    
It also makes it easy to create new time ranges. The functions next() and previous()
do an inplace update to the object by either adding or subtracting the same time interval
. This could be useful if you need to step through a number of time ranges. For example,
if you needed time ranges that spanned 30 minutes over a period of 4 hours you could do: ::

    for a in range(8): print(time_range.next())
    
Check out the code reference for the time range object for more information.