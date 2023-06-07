.. _sunpy-topic-guide-timeseries-metadata:

A Detailed Look at the TimeSeries Metadata
==========================================

TimeSeries store metadata in a `~sunpy.timeseries.TimeSeriesMetaData` object, this object is designed to be able to store multiple basic `~sunpy.util.metadata.MetaDict` (case-insensitive ordered dictionary) objects and able to identify the relevant metadata for a given cell in the data.
This enables a single TimeSeries to be created by combining/concatenating multiple TimeSeries source files together into one and to keep a reliable track of all the metadata relevant to each cell, column or row.
The metadata can be accessed by:

.. code-block:: python

    >>> import sunpy.timeseries
    >>> import sunpy.data.sample  # doctest: +REMOTE_DATA

    >>> my_timeseries = sunpy.timeseries.TimeSeries(sunpy.data.sample.GOES_XRS_TIMESERIES)  # doctest: +REMOTE_DATA
    >>> meta = my_timeseries.meta # doctest: +REMOTE_DATA

You can easily get an overview of the metadata, this will show you a basic representation of the metadata entries that are relevant to this TimeSeries.

.. code-block:: python

    >>> meta # doctest: +REMOTE_DATA
    |-------------------------------------------------------------------------------------------------|
    |TimeRange                  | Columns         | Meta                                              |
    |-------------------------------------------------------------------------------------------------|
    |2011-06-06T23:59:59.961999 | xrsa            | simple: True                                      |
    |            to             | xrsb            | bitpix: 8                                         |
    |2011-06-07T23:59:57.631999 |                 | naxis: 0                                          |
    |                           |                 | extend: True                                      |
    |                           |                 | date: 26/06/2012                                  |
    |                           |                 | numext: 3                                         |
    |                           |                 | telescop: GOES 15                                 |
    |                           |                 | instrume: X-ray Detector                          |
    |                           |                 | object: Sun                                       |
    |                           |                 | origin: SDAC/GSFC                                 |
    |                           |                 | ...                                               |
    |-------------------------------------------------------------------------------------------------|
    <BLANKLINE>

The data within a `~sunpy.timeseries.TimeSeriesMetaData` object is stored as a list of tuples, each tuple representing the metadata from a source file or provided when creating the time series.
The tuple will contain a `~sunpy.time.TimeRange` telling us which rows the metadata applies to, a list of column name strings for which the metadata applies to and finally a `~sunpy.util.metadata.MetaDict` object for storing the key/value pairs of the metadata itself.
Each time a TimeSeries is concatenated to the original a new set of rows and/or columns will be added to the `~pandas.DataFrame` and a new entry will be added into the metadata.
Note that entries are ordered chronologically based on `~sunpy.time.timerange.TimeRange.start` and generally it's expected that no two TimeSeries will overlap on both columns and time range.
For example, it is not good practice for alternate row values in a single column to be relevant to different metadata entries as this would make it impossible to uniquely identify the metadata relevant to each cell.

If you want the string that is printed then you can use the `~sunpy.timeseries.TimeSeriesMetaData.to_string` method.
This has the advantage of having optional keyword arguments that allows you to set the depth (number of rows for each entry) and width (total number of characters wide) to better fit your output.
For example:

.. code-block:: python

    >>> meta.to_string(depth=20, width=99) # doctest: +REMOTE_DATA
    "|-------------------------------------------------------------------------------------------------|\n|TimeRange                  | Columns         | Meta                                              |\n|-------------------------------------------------------------------------------------------------|\n|2011-06-06T23:59:59.961999 | xrsa            | simple: True                                      |\n|            to             | xrsb            | bitpix: 8                                         |\n|2011-06-07T23:59:57.631999 |                 | naxis: 0                                          |\n|                           |                 | extend: True                                      |\n|                           |                 | date: 26/06/2012                                  |\n|                           |                 | numext: 3                                         |\n|                           |                 | telescop: GOES 15                                 |\n|                           |                 | instrume: X-ray Detector                          |\n|                           |                 | object: Sun                                       |\n|                           |                 | origin: SDAC/GSFC                                 |\n|                           |                 | date-obs: 07/06/2011                              |\n|                           |                 | time-obs: 00:00:00.000                            |\n|                           |                 | date-end: 07/06/2011                              |\n|                           |                 | time-end: 23:59:57.632                            |\n|                           |                 | comment: Energy band information given in extensio|\n|                           |                 | history:                                          |\n|                           |                 | keycomments: {'SIMPLE': 'Written by IDL:  Tue Jun |\n|-------------------------------------------------------------------------------------------------|\n"


Similar to the TimeSeries, the metadata has some properties for convenient access to the global metadata details, including `~sunpy.timeseries.TimeSeriesMetaData.columns` as a list of strings,  and `~sunpy.timeseries.TimeSeriesMetaData.time_range` of the data.
Beyond this, there are properties to get lists of details for all the entries in the `~sunpy.timeseries.TimeSeriesMetaData` object, including `~sunpy.timeseries.TimeSeriesMetaData.timeranges`, `~sunpy.timeseries.TimeSeriesMetaData.columns` (as a list of string column names) and `~sunpy.timeseries.TimeSeriesMetaData.metas`.
Similar to TimeSeries objects you can `~sunpy.timeseries.TimeSeriesMetaData.concatenate` `~sunpy.timeseries.TimeSeriesMetaData` objects, but generally you won't need to do this as it is done automatically when actioned on the TimeSeries.
Note that when truncating a `~sunpy.timeseries.TimeSeriesMetaData` object you will remove any entries outside of the given `~sunpy.time.TimeRange`.
You can also `~sunpy.timeseries.TimeSeriesMetaData.append` a new entry (as a tuple or list), which will add the entry in the correct chronological position.
It is frequently necessary to locate the metadata for a given column, row or cell which can be uniquely identified by both, to do this you can use the `~sunpy.timeseries.TimeSeriesMetaData.find` method, by adding "colname" and/or time/row keyword arguments you get a `~sunpy.timeseries.TimeSeriesMetaData` object returned which contains only the relevant entries.
You can then use the `~sunpy.timeseries.TimeSeriesMetaData.metas` property to get a list of just the relevant `~sunpy.util.metadata.MetaDict` objects.
For example:

.. code-block:: python

    >>> tsmd_return = my_timeseries.meta.find(colname='xrsa', time='2012-06-01 00:00:33.904999') # doctest: +REMOTE_DATA
    >>> tsmd_return.metas # doctest: +REMOTE_DATA
    []

Note, the ``colname`` and ``time`` keywords are optional, but omitting both just returns the original.
A common use case for the metadata is to find out the instrument(s) that gathered the data and in this case you can use the `~sunpy.timeseries.TimeSeriesMetaData.get` method.
This method takes a single key string or list of key strings with the optional filters and will search for any matching values.
This method returns another `~sunpy.timeseries.TimeSeriesMetaData` object, but removes all unwanted key/value pairs.
The result can be converted into a simple list of strings using the `~sunpy.timeseries.TimeSeriesMetaData.values` method:

.. code-block:: python

    >>> tsmd_return = my_timeseries.meta.get('telescop', colname='xrsa') # doctest: +REMOTE_DATA
    >>> tsmd_return.values() # doctest: +REMOTE_DATA
    ['GOES 15']

Note `~sunpy.timeseries.TimeSeriesMetaData.values` removes duplicate strings and sorts the returned list.
You can update the values for these entries efficiently using the `~sunpy.timeseries.TimeSeriesMetaData.update` method which takes a dictionary argument and updates the values to each of the dictionaries that match the given ``colname`` and ``time`` filters, for example:

.. code-block:: python

    >>> my_timeseries.meta.update({'telescop': 'G15'}, colname='xrsa', overwrite=True) # doctest: +REMOTE_DATA

Here we have to specify the ``overwrite=True`` keyword parameter to allow us to overwrite values for keys already present in the `~sunpy.util.metadata.MetaDict` objects, this helps protect the integrity of the original metadata and without this set (or with it set to `False`) you can still add new key/value pairs.
Note that the `~sunpy.util.metadata.MetaDict` objects are both case-insensitive for key strings and have ordered entries, where possible the order is preserved when updating values.
