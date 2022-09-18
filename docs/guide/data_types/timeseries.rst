.. _timeseries_guide:

**********************
Time series data guide
**********************

Time series data are a fundamental part of many data analysis projects in heliophysics as well as other areas.
**sunpy** provides a TimeSeries object to handle this type of data.
Much like the `~sunpy.map.Map` object, `~sunpy.timeseries.TimeSeries` can load generic filetypes, and recognizes data from nine specific sources to provide instrument specific data loading and plotting capabilities.
For more information about TimeSeries, and what file types and data sources are supported, check out :doc:`/code_ref/timeseries`.

:ref:`creating-timeseries` describes how to create a TimeSeries object from single or multiple observational sources.
:ref:`inspecting-timeseries` describes how to examine the data and metadata.
:ref:`plotting-timeseries` outlines the basics of how **sunpy** will plot a TimeSeries object.
:ref:`manipulating-timeseries` describes how to modify data, truncate a TimeSeries, down- and up-sample data, concatenate data, and create an Astropy Table from a TimeSeries.
:ref:`custom-timeseries` details how to create a TimeSeries object from a Pandas DataFrame or an Astropy Table.
Lastly, :ref:`ts-metadata` describes how to view and extract information from the metadata.

.. _creating-timeseries:

1. Creating a TimeSeries
========================

A TimeSeries object can be created from local files.
For convenience, **sunpy** can download several example timeseries of observational data.
These files have names like ``sunpy.data.sample.EVE_TIMESERIES`` and ``sunpy.data.sample.GOES_XRS_TIMESERIES``.
To create the sample `sunpy.timeseries.sources.goes.XRSTimeSeries`, type the following into your interactive Python shell:

.. code-block:: python

    >>> import sunpy.timeseries as ts
    >>> import sunpy.data.sample  # doctest: +REMOTE_DATA
    >>> my_timeseries = ts.TimeSeries(sunpy.data.sample.GOES_XRS_TIMESERIES)  # doctest: +REMOTE_DATA

This is calling the `~sunpy.timeseries.TimeSeries` factory to create a time series from a GOES XRS FITS file.

The variable ``my_timeseries`` is a `~sunpy.timeseries.GenericTimeSeries` object.
To create one from a local GOES/XRS FITS file try the following:

.. code-block:: python

    >>> my_timeseries = ts.TimeSeries('/mydirectory/myts.fits', source='XRS')   # doctest: +SKIP

**sunpy** will attempt to detect automatically the instrument source for most FITS files.
However timeseries data are stored in a variety of file types (FITS, txt, csv, CDF), and so it is not always possible to detect the source.
**sunpy** ships with a number of known instrumental sources, and can also load CDF files that conform to the `Space Physics Guidelines for CDF <https://spdf.gsfc.nasa.gov/sp_use_of_cdf.html>`__.
If you would like **sunpy** to include another instrumental source see the `Newcomers' Guide <https://docs.sunpy.org/en/latest/dev_guide/contents/newcomers.html>`__.

The `~sunpy.timeseries.TimeSeries` factory has the ability to create a list of TimeSeries objects using a list of filepaths, a folder or a glob, for example:

.. code-block:: python

    >>> my_ts_list = ts.TimeSeries('filepath1', 'filepath2', source='XRS')   # doctest: +SKIP
    >>> my_ts_list = ts.TimeSeries('/goesdirectory/', source='XRS')   # doctest: +SKIP
    >>> my_ts_list = ts.TimeSeries(glob, source='XRS')   # doctest: +SKIP

When manually specifing the source this functionality will only work with files from the same single source, generating a source specific child of the `~sunpy.timeseries.GenericTimeSeries` class such as the `~sunpy.timeseries.sources.goes.XRSTimeSeries` above.

Instead of creating a list of one TimeSeries object per file, you can create a single time series from multiple files using the keyword argument ``concatenate=True``:

.. code-block:: python

    >>> my_ts = ts.TimeSeries(sunpy.data.sample.GOES_XRS_TIMESERIES, sunpy.data.sample.GOES_XRS_TIMESERIES, source='XRS', concatenate=True)  # doctest: +REMOTE_DATA

Again these must all be from the same source if the ``source`` keyword is explicitly specified.
The `.GenericTimeSeries.concatenate` method can be used to make a single time series from multiple TimeSeries from different sources if they are already in the form of TimeSeries objects.

.. _inspecting-timeseries:

2. Inspecting TimeSeries & Accessing the Data
=============================================

A TimeSeries object holds both data as well as metadata and unit data.
For a quick look at a TimeSeries, type:

.. code-block:: python

    >>> my_timeseries  # doctest: +REMOTE_DATA
    <sunpy.timeseries.sources.goes.XRSTimeSeries object at ...>
    SunPy TimeSeries
    ----------------
    Observatory:		 GOES-15
    Instrument:		 <a href=https://www.swpc.noaa.gov/products/goes-x-ray-flux target="_blank">X-ray Detector</a>
    Channel(s):		 xrsa<br>xrsb
    Start Date:		 2011-06-07 00:00:00
    End Date:		 2011-06-07 23:59:58
    Center Date:		 2011-06-07 11:59:58
    Resolution:		 2.048 s
    Samples per Channel:		 42177
    Data Range(s):		 xrsa   3.64E-06<br>xrsb   2.54E-05
    Units:		 W / m2
                                           xrsa          xrsb
    2011-06-06 23:59:59.961999893  1.000000e-09  1.887100e-07
    2011-06-07 00:00:02.008999944  1.000000e-09  1.834600e-07
    2011-06-07 00:00:04.058999896  1.000000e-09  1.860900e-07
    2011-06-07 00:00:06.104999900  1.000000e-09  1.808400e-07
    2011-06-07 00:00:08.151999950  1.000000e-09  1.860900e-07
    ...                                     ...           ...
    2011-06-07 23:59:49.441999912  1.000000e-09  1.624800e-07
    2011-06-07 23:59:51.488999844  1.000000e-09  1.624800e-07
    2011-06-07 23:59:53.538999915  1.000000e-09  1.598500e-07
    2011-06-07 23:59:55.584999919  1.000000e-09  1.624800e-07
    2011-06-07 23:59:57.631999850  1.000000e-09  1.598500e-07
    <BLANKLINE>
    [42177 rows x 2 columns]

This shows a table of information taken from the metadata and a preview of your data.
If you execute this command in a Jupyter Notebook, a rich HTML version of this quick look will be shown that includes plots of the data.
Alternatively, the :func:`~sunpy.timeseries.GenericTimeSeries.quicklook` command will show the HTML view in your default browser.
The metadata for the time series is accessed by:

.. code-block:: python

    >>> my_timeseries.meta # doctest: +REMOTE_DATA
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

This references the `~sunpy.timeseries.TimeSeriesMetaData` object with the header information as read from the source files.
A word of caution: many data sources provide little to no meta data so this variable might be empty.
The meta data is described in more detail later in this guide.
Similarly there are properties for getting `~sunpy.timeseries.GenericTimeSeries.columns` as a list of strings, `~sunpy.timeseries.GenericTimeSeries.time` values and `~sunpy.timeseries.GenericTimeSeries.time_range` of the data.

To get a column of the data use the `~sunpy.timeseries.GenericTimeSeries.quantity` method:

.. code-block:: python

    >>> my_timeseries.quantity('xrsa') # doctest: +REMOTE_DATA
    <Quantity [1.e-09, 1.e-09, 1.e-09, ..., 1.e-09, 1.e-09, 1.e-09] W / m2>

.. _plotting-timeseries:

3. Plotting TimeSeries
======================

The **sunpy** TimeSeries object has its own built-in plot methods so that it is easy to quickly view your time series.
To create a plot just type:

.. plot::
    :include-source:

    import sunpy.timeseries as ts
    import sunpy.data.sample

    ts = ts.TimeSeries(sunpy.data.sample.GOES_XRS_TIMESERIES, source='XRS')
    ts.peek()

This will open a Matplotlib plot on your screen.
If you want to save this to a PNG file you can do so from the Matplotlib GUI.

In addition, to enable users to modify the plot it is possible to use the `~sunpy.timeseries.GenericTimeSeries.plot` command.
This makes it possible to use the **sunpy** plot as the foundation for a more complicated figure:

.. plot::
   :include-source:

   import matplotlib.pyplot as plt

   import sunpy.timeseries as ts
   import sunpy.data.sample

   ts = ts.TimeSeries(sunpy.data.sample.GOES_XRS_TIMESERIES, source='XRS')
   fig, ax = plt.subplots()
   ts.plot(axes=ax)
   # Add code to modify the figure here if desired
   fig.savefig('figure.png')

.. _manipulating-timeseries:

4. Manipulating TimeSeries
==========================

4.1 Modifying the Data
----------------------
TimeSeries provides the `~sunpy.timeseries.GenericTimeSeries.add_column` method which will either add a new column or update a current column if the colname is already present.
This can take numpy array or preferably an Astropy `~astropy.units.quantity.Quantity` value.
For example:

.. code-block:: python

    >>> values = my_timeseries.quantity('xrsa') * 2
    >>> my_timeseries = my_timeseries.add_column('xrsa*2', values)
    >>> my_timeseries.columns
    ['xrsa', 'xrsb', 'xrsa*2']

Adding a column is not done in place, but instead returns a new TimeSeries with the new column added.
Note that the values will be converted into the column units if an Astropy `~astropy.units.quantity.Quantity` is given.
Caution should be taken when adding a new column because this column won't have any associated MetaData entry.

4.2 Truncating a TimeSeries
---------------------------

It is often useful to truncate an existing TimeSeries object to retain a specific time range.
This is easily achieved by using the `~sunpy.timeseries.GenericTimeSeries.truncate` method.
For example, to trim our GOES data into a period of interest use:

.. code-block:: python

    >>> from sunpy.time import TimeRange
    >>> tr = TimeRange('2012-06-01 05:00', '2012-06-01 06:30')
    >>> my_timeseries_trunc = my_timeseries.truncate(tr)

This takes a number of different arguments, such as the start and end dates (as datetime or string objects) or a `~sunpy.time.TimeRange` as used above.
Note that the truncated TimeSeries will have a truncated `~sunpy.timeseries.TimeSeriesMetaData` object, which may include dropping metadata entries for data totally cut out from the TimeSeries.
If you want to truncate using slice-like values you can, for example taking every 2nd value from 0 to 10000 can be done using:

.. code-block:: python

    >>> my_timeseries_trunc = my_timeseries.truncate(0, 100000, 2)

4.3 More complicated timeseries operations
------------------------------------------
If you want to do any more complicated analysis on a TimeSeries, we recommend converting it to a `pandas.DataFrame` object first.
Although this conversion will use the unit information and metadata, pandas has a wide array of methods that can be used e.g. for resampling data.
As an example to downsample you can do:

.. code-block:: python

    >>> downsampled_dataframe = my_timeseries_trunc.to_dataframe().resample('10T').mean()

Here ``10T`` means sample every 10 minutes and 'mean' is the method used to combine the data in each 10 minute bin.
See the `pandas` documentation for more details on other functionality they offer for timeseries analysis.

4.4 Concatenating TimeSeries
----------------------------
It's common to want to combine a number of TimeSeries together into a single TimeSeries.
In the simplest scenario this is to combine data from a single source over several time ranges, for example if you wanted to combine the daily GOES data to get a week or more of constant data in one TimeSeries.
This can be performed using the TimeSeries factory with the ``concatenate=True`` keyword argument:

.. code-block:: python

    >>> concatenated_timeseries = sunpy.timeseries.TimeSeries(filepath1, filepath2, source='XRS', concatenate=True)  # doctest: +SKIP

Note, you can list any number of files, or a folder or use a glob to select the input files to be concatenated.
It is possible to concatenate two TimeSeries after creating them with the factory using the `~sunpy.timeseries.GenericTimeSeries.concatenate` method.
For example:

.. code-block:: python

    >>> concatenated_timeseries = goes_timeseries_1.concatenate(goes_timeseries_2) # doctest: +SKIP

This will result in a TimeSeries identical to if you used the factory to create it in one step.
A limitation of the TimeSeries class is that often it is not easy to determine the source observatory/instrument of a file, generally because the file formats used vary depending on the scientific working groups, thus some sources need to be explicitly stated (as a keyword argument) and so it is not possible to concatenate files from multiple sources with the factory.
To do this you can still use the `~sunpy.timeseries.GenericTimeSeries.concatenate` method, which will create a new TimeSeries with all the rows and columns of the source and concatenated TimeSeries in one:

.. code-block:: python

    >>> eve_ts = ts.TimeSeries(sunpy.data.sample.EVE_TIMESERIES, source='eve')
    >>> goes_ts = ts.TimeSeries(sunpy.data.sample.GOES_XRS_TIMESERIES)
    >>> concatenated_timeseries = goes_ts.concatenate(eve_ts)

Note that the more complex `~sunpy.timeseries.TimeSeriesMetaData` object now has 2 entries and shows details on both:

.. code-block:: python

    >>> concatenated_timeseries.meta
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
    |2011-06-07T00:00:00.000000 | XRS-B proxy     | data_list: 20110607_EVE_L0CS_DIODES_1m.txt        |
    |            to             | XRS-A proxy     | created: Tue Jun  7 23:59:10 2011 UTC             |
    |2011-06-07T23:59:00.000000 | SEM proxy       | origin: SDO/EVE Science Processing and Operations |
    |                           | 0.1-7ESPquad    | units: W/m^2 for irradiance, dark is counts/(0.25s|
    |                           | 17.1ESP         | source: SDO-EVE ESP and MEGS-P instruments, http:/|
    |                           | 25.7ESP         | product: Level 0CS, 1-minute averaged SDO-EVE Sola|
    |                           | 30.4ESP         | version: 2.1, code updated 2011-May-12            |
    |                           | 36.6ESP         | missing data: -1.00e+00                           |
    |                           | darkESP         | hhmm: hour and minute in UT                       |
    |                           | 121.6MEGS-P     | xrs-b proxy: a model of the expected XRS-B 0.1-0.8|
    |                           | ...             | ...                                               |
    |-------------------------------------------------------------------------------------------------|
    <BLANKLINE>

The metadata object is described in more detail in the next section.


4.5 Creating an Astropy Table from a TimeSeries
-----------------------------------------------

If you want to take the data from your TimeSeries and use it as a `~astropy.table.Table` this can be done using the `~sunpy.timeseries.GenericTimeSeries.to_table` method.
For example:

.. code-block:: python

    >>> table = my_timeseries_trunc.to_table()

Note that this `~astropy.table.Table` will contain a mixin column for containing the Astropy `~astropy.time.Time` object representing the index, it will also add the relevant units to the columns.
One of the most useful reasons for doing this is that Astropy `~sunpy.timeseries.GenericTimeSeries.to_table` objects have some very nice options for viewing the data, including the basic console view:

.. code-block:: python

    >>> table
    <Table length=21089>
                 date               xrsa     xrsb     xrsa*2
                                   W / m2   W / m2    W / m2
            datetime64[ns]        float32  float32   float32
    ----------------------------- ------- ---------- -------
    2011-06-06T23:59:59.961999893   1e-09 1.8871e-07   2e-09
    2011-06-07T00:00:04.058999896   1e-09 1.8609e-07   2e-09
    2011-06-07T00:00:08.151999950   1e-09 1.8609e-07   2e-09
    2011-06-07T00:00:12.248999953   1e-09 1.8609e-07   2e-09
    2011-06-07T00:00:16.344999909   1e-09 1.8084e-07   2e-09
    2011-06-07T00:00:20.441999912   1e-09 1.8084e-07   2e-09
    2011-06-07T00:00:24.534999847   1e-09 1.8084e-07   2e-09
    2011-06-07T00:00:28.631999850   1e-09 1.8346e-07   2e-09
    2011-06-07T00:00:32.728999853   1e-09 1.8346e-07   2e-09
                              ...     ...        ...     ...
    2011-06-07T23:59:20.768999934   1e-09  1.651e-07   2e-09
    2011-06-07T23:59:24.864999890   1e-09 1.5985e-07   2e-09
    2011-06-07T23:59:28.961999893   1e-09 1.5985e-07   2e-09
    2011-06-07T23:59:33.058999896   1e-09 1.6248e-07   2e-09
    2011-06-07T23:59:37.151999950   1e-09 1.6248e-07   2e-09
    2011-06-07T23:59:41.248999953   1e-09 1.5985e-07   2e-09
    2011-06-07T23:59:45.344999909   1e-09 1.5723e-07   2e-09
    2011-06-07T23:59:49.441999912   1e-09 1.6248e-07   2e-09
    2011-06-07T23:59:53.538999915   1e-09 1.5985e-07   2e-09
    2011-06-07T23:59:57.631999850   1e-09 1.5985e-07   2e-09

and the more sophisticated browser view using the `~astropy.table.Table.show_in_browser` method:

.. code-block:: python

    >>> table.show_in_browser(jsviewer=True)  # doctest: +SKIP

For further details about editing Astropy tables you can read the `Astropy documentation website <https://docs.astropy.org/en/stable/table/>`_.

.. _custom-timeseries:

5. Creating Custom TimeSeries
=============================

Sometimes you will have data that you want to create into a TimeSeries.
You can use the factory to create a `~sunpy.timeseries.GenericTimeSeries` from a variety of data sources currently including `pandas.DataFrame` and `astropy.table.Table`.

5.1 Creating a TimeSeries from a Pandas DataFrame
-------------------------------------------------

A TimeSeries object must be supplied with some data when it is created.
The data can either be in your current Python session, in a local file, or in a remote file.
Let's create some data and pass it into a TimeSeries object:

.. code-block:: python

    >>> import numpy as np
    >>> intensity = np.sin(np.arange(0, 12 * np.pi, ((12 * np.pi) / (24*60))))

This creates a basic numpy array of values representing a sine wave.
We can use this array along with a suitable time storing object (such as Astropy `~astropy.time` or a list of `datetime` objects) to make a Pandas `~pandas.DataFrame`.
A suitable list of times must contain the same number of values as the data, this can be created using:

.. code-block:: python

    >>> import datetime
    >>> base = datetime.datetime.today()
    >>> times = [base - datetime.timedelta(minutes=x) for x in range(24*60, 0, -1)]

The Pandas `~pandas.DataFrame` will use the dates list as the index:

.. code-block:: python

    >>> from pandas import DataFrame
    >>> data = DataFrame(intensity, index=times, columns=['intensity'])

This `~pandas.DataFrame` can then be used to construct a TimeSeries:

.. code-block:: python

    >>> import sunpy.timeseries as ts
    >>> import astropy.units as u
    >>> header = {'key': 'value'}
    >>> units = {'intensity': u.W/u.m**2}
    >>> ts_custom = ts.TimeSeries(data, header, units)

5.2 Creating Custom TimeSeries from an Astropy Table
----------------------------------------------------

A Pandas `~pandas.DataFrame` is the underlying object used to store the data within a TimeSeries, so the above example is the most lightweight to create a custom TimeSeries, but being scientific data it will often be more convenient to use an Astropy `~astropy.table.Table` and let the factory convert this.
An advantage of this method is it allows you to include metadata and Astropy `~astropy.units.quantity.Quantity` values, which are both supported in tables, without additional arguments.
For example:

.. code-block:: python

    >>> import datetime
    >>> from astropy.time import Time
    >>> import astropy.units as u
    >>> from astropy.table import Table

    >>> base = datetime.datetime.today()
    >>> times = [base - datetime.timedelta(minutes=x) for x in range(24*60, 0, -1)]
    >>> intensity = u.Quantity(np.sin(np.arange(0, 12 * np.pi, ((12 * np.pi) / (24*60)))), u.W/u.m**2)
    >>> tbl_meta = {'t_key':'t_value'}
    >>> table = Table([times, intensity], names=['time', 'intensity'], meta=tbl_meta)
    >>> table.add_index('time')
    >>> ts_table = ts.TimeSeries(table)

Note that due to the properties of the `~astropy.time.Time` object, this will be a mixin column which since it is a single object, limits the versatility of the `~astropy.table.Table` a little.
For more on mixin columns see the `Astropy docs <https://docs.astropy.org/en/stable/table/mixin_columns.html>`__.
The units will be taken from the table quantities for each column, the metadata will simply be the table.meta dictionary.
You can also explicitly add metadata and units, these will be added to the relevant dictionaries using the dictionary update method, with the explicit user-given values taking precedence:

.. code-block:: python

    >>> from sunpy.util.metadata import MetaDict
    >>> from collections import OrderedDict
    >>> import astropy.units as u

    >>> meta = MetaDict({'key':'value'})
    >>> units = OrderedDict([('intensity', u.W/u.m**2)])
    >>> ts_table = ts.TimeSeries(table, meta, units)

.. _ts-metadata:

6. A Detailed Look at the Metadata
==================================

TimeSeries store metadata in a `~sunpy.timeseries.TimeSeriesMetaData` object, this object is designed to be able to store multiple basic `~sunpy.util.metadata.MetaDict` (case-insensitive ordered dictionary) objects and able to identify the relevant metadata for a given cell in the data.
This enables a single TimeSeries to be created by combining/concatenating multiple TimeSeries source files together into one and to keep a reliable track of all the metadata relevant to each cell, column or row.
The metadata can be accessed by:

.. code-block:: python

    >>> meta = my_timeseries.meta

You can easily get an overview of the metadata, this will show you a basic representation of the metadata entries that are relevant to this TimeSeries.

.. code-block:: python

    >>> meta
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

The data within a `~sunpy.timeseries.TimeSeriesMetaData` object is stored as a list of tuples, each tuple representing the metadata from a source file or timeseries.
The tuple will contain a `~sunpy.time.TimeRange` telling us which rows the metadata applies to, a list of column name strings for which the metadata applies to and finally a `~sunpy.util.metadata.MetaDict` object for storing the key/value pairs of the metadata itself.
Each time a TimeSeries is concatenated to the original a new set of rows and/or columns will be added to the `~pandas.DataFrame` and a new entry will be added into the metadata.
Note that entries are ordered chronologically based on
`~sunpy.time.timerange.TimeRange.start` and generally it's expected that no two
TimeSeries will overlap on both columns and time range.
For example, it is not good practice for alternate row values in a single column to be relevant to different metadata entries as this would make it impossible to uniquely identify the metadata relevant to each cell.

If you want the string that's printed then you can use the `~sunpy.timeseries.TimeSeriesMetaData.to_string` method.
This has the advantage of having optional keyword arguments that allows you to set the depth (number of rows for each entry) and width (total number of characters wide) to better fit your output.
For example:

.. code-block:: python

    >>> meta_str = meta.to_string(depth = 20, width=99)

Similar to the TimeSeries, the metadata has some properties for convenient access to the global metadata details, including `~sunpy.timeseries.TimeSeriesMetaData.columns` as a list of strings,  and `~sunpy.timeseries.TimeSeriesMetaData.time_range` of the data.
Beyond this, there are properties to get lists of details for all the entries in the `~sunpy.timeseries.TimeSeriesMetaData` object, including `~sunpy.timeseries.TimeSeriesMetaData.timeranges`, `~sunpy.timeseries.TimeSeriesMetaData.columns` (as a list of string column names) and `~sunpy.timeseries.TimeSeriesMetaData.metas`.
Similar to TimeSeries objects you can `~sunpy.timeseries.TimeSeriesMetaData.concatenate` `~sunpy.timeseries.TimeSeriesMetaData` objects, but generally you won't need to do this as it is done automatically when actioned on the TimeSeries.
Note that when truncating a `~sunpy.timeseries.TimeSeriesMetaData` object you will remove any entries outside of the given `~sunpy.time.TimeRange`.
You can also `~sunpy.timeseries.TimeSeriesMetaData.append` a new entry (as a tuple or list), which will add the entry in the correct chronological position.
It is frequently necessary to locate the metadata for a given column, row or cell which can be uniquely identified by both, to do this you can use the `~sunpy.timeseries.TimeSeriesMetaData.find` method, by adding colname and/or time/row keyword arguments you get a `~sunpy.timeseries.TimeSeriesMetaData` object returned which contains only the relevant entries.
You can then use the `~sunpy.timeseries.TimeSeriesMetaData.metas` property to get a list of just the relevant `~sunpy.util.metadata.MetaDict` objects.
For example:

.. code-block:: python

    >>> tsmd_return = my_timeseries.meta.find(colname='xrsa', time='2012-06-01 00:00:33.904999')
    >>> tsmd_return.metas
    []

Note, the colname and time filters are optional, but omitting both filters just returns an identical `~sunpy.timeseries.TimeSeriesMetaData` object to the TimeSeries original.
A common use case for the metadata is to find out the instrument/s that gathered the data and in this case you can use the `~sunpy.timeseries.TimeSeriesMetaData.get` method.
This method takes a single key string or list of key strings with the optional filters and will search for any matching values.
This method returns another `~sunpy.timeseries.TimeSeriesMetaData` object, but removes all unwanted key/value pairs.
The result can be converted into a simple list of strings using the `~sunpy.timeseries.TimeSeriesMetaData.values` method:

.. code-block:: python

    >>> tsmd_return = my_timeseries.meta.get('telescop', colname='xrsa')
    >>> tsmd_return.values()
    ['GOES 15']

Note `~sunpy.timeseries.TimeSeriesMetaData.values` removes duplicate strings and sorts the returned list.
You can update the values for these entries efficiently using the `~sunpy.timeseries.TimeSeriesMetaData.update` method which takes a dictionary argument and updates the values to each of the dictionaries that match the given colname and time filters, for example:

.. code-block:: python

    >>> my_timeseries.meta.update({'telescop': 'G15'}, colname='xrsa', overwrite=True)

Here we have to specify the ``overwrite=False`` keyword parameter to allow us to overwrite values for keys already present in the `~sunpy.util.metadata.MetaDict` objects, this helps protect the integrity of the original metadata and without this set (or with it set to False) you can still add new key/value pairs.
Note that the `~sunpy.util.metadata.MetaDict` objects are both case-insensitive for key strings and have ordered entries, where possible the order is preserved when updating values.
