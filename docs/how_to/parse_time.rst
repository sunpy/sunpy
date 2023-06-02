.. _sunpy-how-to-parse-times-with-parse-time:

*****************************************
Parse times with `~sunpy.time.parse_time`
*****************************************

.. code-block:: python

    >>> import time
    >>> from datetime import date, datetime

    >>> import numpy as np
    >>> import pandas

    >>> from sunpy.time import parse_time

The following examples show how to use `sunpy.time.parse_time` to parse various time formats, including both strings and objects, into an `astropy.time.Time` object.

Strings
=======

.. code-block:: python

    >>> parse_time('1995-12-31 23:59:60')
    <Time object: scale='utc' format='isot' value=1995-12-31T23:59:60.000>

This also works with the ``scale=`` keyword argument (See `this list: <https://docs.astropy.org/en/stable/time/#time-scale>`__ for the list of all allowed scales):

.. code-block:: python

    >>> parse_time('2012:124:21:08:12', scale='tai')
    <Time object: scale='tai' format='isot' value=2012-05-03T21:08:12.000>

Tuples
======

.. code-block:: python

    >>> parse_time((1998, 11, 14))
    <Time object: scale='utc' format='isot' value=1998-11-14T00:00:00.000>
    >>> parse_time((2001, 1, 1, 12, 12, 12, 8899))
    <Time object: scale='utc' format='isot' value=2001-01-01T12:12:12.009>

`time.struct_time`
==================

.. code-block:: python

    >>> parse_time(time.gmtime(0))
    <Time object: scale='utc' format='isot' value=1970-01-01T00:00:00.000>

`datetime.datetime` and `datetime.date`
=======================================

.. code-block:: python

    >>> parse_time(datetime(1990, 10, 15, 14, 30))
    <Time object: scale='utc' format='datetime' value=1990-10-15 14:30:00>
    >>> parse_time(date(2023, 4, 22))
    <Time object: scale='utc' format='iso' value=2023-04-22 00:00:00.000>

`pandas` time objects
=====================

`pandas.Timestamp`, `pandas.Series` and `pandas.DatetimeIndex`

.. code-block:: python

    >>> parse_time(pandas.Timestamp(datetime(1966, 2, 3)))
    <Time object: scale='utc' format='datetime64' value=1966-02-03T00:00:00.000000000>
    >>> parse_time(pandas.Series([[datetime(2012, 1, 1, 0, 0), datetime(2012, 1, 2, 0, 0)],
    ...                           [datetime(2012, 1, 3, 0, 0), datetime(2012, 1, 4, 0, 0)]]))
    <Time object: scale='utc' format='datetime' value=[[datetime.datetime(2012, 1, 1, 0, 0) datetime.datetime(2012, 1, 2, 0, 0)]
                                                       [datetime.datetime(2012, 1, 3, 0, 0) datetime.datetime(2012, 1, 4, 0, 0)]]>
    >>> parse_time(pandas.DatetimeIndex([datetime(2012, 1, 1, 0, 0),
    ...                                  datetime(2012, 1, 2, 0, 0),
    ...                                  datetime(2012, 1, 3, 0, 0),
    ...                                  datetime(2012, 1, 4, 0, 0)]))
    <Time object: scale='utc' format='datetime' value=[datetime.datetime(2012, 1, 1, 0, 0)
                                                       datetime.datetime(2012, 1, 2, 0, 0)
                                                       datetime.datetime(2012, 1, 3, 0, 0)
                                                       datetime.datetime(2012, 1, 4, 0, 0)]>

`numpy.datetime64`
==================

.. code-block:: python

    >>> parse_time(np.datetime64('2014-02-07T16:47:51.008288123'))
    <Time object: scale='utc' format='isot' value=2014-02-07T16:47:51.008>
    >>> parse_time(np.array(['2014-02-07T16:47:51.008288123', '2014-02-07T18:47:51.008288123'],
    ...                     dtype='datetime64'))
    <Time object: scale='utc' format='isot' value=['2014-02-07T16:47:51.008' '2014-02-07T18:47:51.008']>

Formats handled by `astropy.time.Time`
======================================

`See this list of all the allowed formats. <https://docs.astropy.org/en/stable/time/#time-format>`__

.. code-block:: python

    >>> parse_time(1234.0, format='jd')
    <Time object: scale='utc' format='jd' value=1234.0>
    >>> parse_time('B1950.0', format='byear_str')
    <Time object: scale='tt' format='byear_str' value=B1950.000>

``anytim`` output
=================

Format output by the ``anytim`` routine in SolarSoft (see the documentation for `~sunpy.time.TimeUTime` for more information):

.. code-block:: python

    >>> parse_time(662738003, format='utime')
    <Time object: scale='utc' format='utime' value=662738003.0>

``anytim2tai`` output
=====================

Format output by the ``anytim2tai`` routine in SolarSoft (see the documentation for `~sunpy.time.TimeTaiSeconds` for more information):

.. code-block:: python

    >>> parse_time(1824441848, format='tai_seconds')
    <Time object: scale='tai' format='tai_seconds' value=1824441848.0>
