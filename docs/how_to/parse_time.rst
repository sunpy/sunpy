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

CDF epochs
==========

Times in `Common Data Format <https://cdf.gsfc.nasa.gov/>`__ (CDF) files, such as the in-situ datasets distributed by `CDAWeb <https://cdaweb.gsfc.nasa.gov/>`__, are stored as plain numbers counted from an epoch.
There are three such types, and the ``format`` to pass to `~sunpy.time.parse_time` for each is:

.. list-table::
   :header-rows: 1

   * - CDF type
     - ``format``
     - Stored as
   * - ``CDF_EPOCH``
     - ``'cdf_epoch'``
     - Milliseconds since 0000-01-01 UTC, as a `float`.
   * - ``CDF_EPOCH16``
     - ``'cdf_epoch16'``
     - Seconds since 0000-01-01 UTC, plus picoseconds within that second in the imaginary part, as a `complex`.
   * - ``CDF_TIME_TT2000``
     - ``'cdf_tt2000'``
     - Nanoseconds since 2000-01-01 12:00:00 TT, as an `int`. Being TT-based, this is the only one of the three that can represent a leap second.

This requires `cdflib <https://cdflib.readthedocs.io/>`__ to be installed, which you will already have if you installed sunpy with the ``timeseries`` extra.

.. doctest-requires:: cdflib

    >>> parse_time(63871286400000.0, format='cdf_epoch')
    <Time object: scale='utc' format='cdf_epoch' value=63871286400000.0>
    >>> parse_time(757339269184000000, format='cdf_tt2000')
    <Time object: scale='tt' format='cdf_tt2000' value=7.57339269184e+17>

Note the scale of the result: because ``CDF_TIME_TT2000`` is counted in TT, the times you get back are on the TT scale, which in 2024 runs 69.184 seconds ahead of UTC.
Convert to UTC before comparing against civil times:

.. doctest-requires:: cdflib

    >>> parse_time(757339269184000000, format='cdf_tt2000').utc.isot
    '2024-01-01T00:00:00.000000000'

Reading a CDF file with `~sunpy.timeseries.TimeSeries` already converts its times for you, so you are most likely to need these formats when the numbers have been separated from the file that defined them.
A time column exported to a CSV, returned by a web API, or listed in someone else's notebook carries no record of the epoch it was counted from.
Given the numbers and which of the three types they are, `~sunpy.time.parse_time` turns them back into times without you having to know how each epoch is defined.

For example, the times below arrived as a list of ``CDF_TIME_TT2000`` values, and can be used as times straight away once parsed:

.. doctest-requires:: cdflib

    >>> times = parse_time([757339269184000000, 757339270184000000, 757339271184000000],
    ...                    format='cdf_tt2000')
    >>> times.utc.isot
    array(['2024-01-01T00:00:00.000000000', '2024-01-01T00:00:01.000000000',
           '2024-01-01T00:00:02.000000000'], dtype='<U29')
    >>> (times[1:] - times[:-1]).sec
    array([1., 1.])
