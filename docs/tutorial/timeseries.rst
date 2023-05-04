.. _timeseries_guide:

**********
Timeseries
**********

Time series data are a fundamental part of many data analysis projects in heliophysics as well as other areas.
**sunpy** provides a TimeSeries object to handle this type of data.
Much like the `~sunpy.map.Map` object, the `~sunpy.timeseries.TimeSeries` object can load generic filetypes, and recognizes data from nine specific sources to provide instrument specific data loading and plotting capabilities.
For more information about TimeSeries, and what file types and data sources are supported, check out :ref:`timeseries_code_ref`.

:ref:`creating-timeseries` describes how to create a TimeSeries object from single or multiple observational sources.
:ref:`inspecting-timeseries` describes how to examine the data and metadata.
:ref:`plotting-timeseries` outlines the basics of how **sunpy** will plot a TimeSeries object.
:ref:`manipulating-timeseries` describes how to modify data, truncate a TimeSeries, down- and up-sample data, concatenate data, and create an Astropy Table from a TimeSeries.
:ref:`custom-timeseries` details how to create a TimeSeries object from a Pandas DataFrame or an Astropy Table.
Lastly, :ref:`ts-metadata` describes how to view and extract information from the metadata.

.. note::

    In this section and in :ref:`map_guide`, we will use the sample data included with sunpy.
    These data are primarily useful for demonstration purposes or simple debugging.
    These files have names like ``sunpy.data.sample.EVE_TIMESERIES`` and ``sunpy.data.sample.GOES_XRS_TIMESERIES`` and are automatically downloaded to your computer as you need them.
    Once downloaded, these sample data files will be paths to their location on your computer.

.. _creating-timeseries:

Creating a TimeSeries
=====================

To create a `sunpy.timeseries.TimeSeries`,:

.. code-block:: python

    >>> import sunpy.timeseries
    >>> import sunpy.data.sample  # doctest: +REMOTE_DATA
    >>> sunpy.data.sample.GOES_XRS_TIMESERIES  # doctest: +REMOTE_DATA
    >>> my_timeseries = sunpy.timeseries.TimeSeries(sunpy.data.sample.GOES_XRS_TIMESERIES)  # doctest: +REMOTE_DATA

In many cases, sunpy will automatically detect the type of the file as well as the instrument associated with it.
In this case, we have a FITS file containing an X-ray light curve as observed by the the XRS instrument on the GOES satellite.

.. note::

    Time series data are stored in a variety of file types (e.g. FITS, csv, CDF), and so it is not always possible to detect the source.
    **sunpy** ships with a number of known instrumental sources, and can also load CDF files that conform to the `Space Physics Guidelines for CDF <https://spdf.gsfc.nasa.gov/sp_use_of_cdf.html>`__.

To make sure this has all worked correctly, we can take a quick look at ``my_timeseries``,

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

This should show a table of information taken from the metadata and a preview of your data.
If you are in a Jupyter Notebook, this will show a rich HTML version that includes plots of the data.
Otherwise, you can use the :meth:`~sunpy.timeseries.GenericTimeSeries.quicklook` method to see this quick-look plots,

.. code-block:: python

    >>> my_timeseries.quicklook()  # doctest: +SKIP

.. generate:: html
    :html_border:

    import sunpy.timeseries
    import sunpy.data.sample
    my_timeseries = sunpy.timeseries.TimeSeries(sunpy.data.sample.GOES_XRS_TIMESERIES)
    print(my_timeseries._repr_html_())

.. _timeseries-data:

TimeSeries Data
===============

To pull out the XRS-A data use the :meth:`~sunpy.timeseries.GenericTimeSeries.quantity` method:

.. code-block:: python

    >>> my_timeseries.quantity('xrsa') # doctest: +REMOTE_DATA
    <Quantity [1.e-09, 1.e-09, 1.e-09, ..., 1.e-09, 1.e-09, 1.e-09] W / m2>

Notice that this is a `~astropy.units.Quantity` object which we discussed in :ref:`units-sunpy`.

.. _inspecting-timeseries:

Inspecting TimeSeries Metadata
===============================

A TimeSeries object holds both data as well as metadata and unit data.
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

.. _plotting-timeseries:

Visualizing TimeSeries
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

Adding Columns
==============

TimeSeries provides the `~sunpy.timeseries.GenericTimeSeries.add_column` method which will either add a new column or update a current column if the colname is already present.
This can take numpy array or preferably an Astropy `~astropy.units.quantity.Quantity` value.
For example:

.. code-block:: python

    >>> values = my_timeseries.quantity('xrsa') * 2 # doctest: +REMOTE_DATA
    >>> my_timeseries = my_timeseries.add_column('xrsa*2', values) # doctest: +REMOTE_DATA
    >>> my_timeseries.columns # doctest: +REMOTE_DATA
    ['xrsa', 'xrsb', 'xrsa*2']

Adding a column is not done in place, but instead returns a new TimeSeries with the new column added.
Note that the values will be converted into the column units if an Astropy `~astropy.units.quantity.Quantity` is given.
Caution should be taken when adding a new column because this column won't have any associated MetaData entry.

Truncating a TimeSeries
=======================

It is often useful to truncate an existing TimeSeries object to retain a specific time range.
This is easily achieved by using the `~sunpy.timeseries.GenericTimeSeries.truncate` method.
For example, to trim our GOES data into a period of interest use:

.. code-block:: python

    >>> from sunpy.time import TimeRange
    >>> tr = TimeRange('2012-06-01 05:00', '2012-06-01 06:30')
    >>> my_timeseries_trunc = my_timeseries.truncate(tr) # doctest: +REMOTE_DATA

This takes a number of different arguments, such as the start and end dates (as datetime or string objects) or a `~sunpy.time.TimeRange` as used above.
Note that the truncated TimeSeries will have a truncated `~sunpy.timeseries.TimeSeriesMetaData` object, which may include dropping metadata entries for data totally cut out from the TimeSeries.
If you want to truncate using slice-like values you can, for example taking every 2nd value from 0 to 10000 can be done using:

.. code-block:: python

    >>> my_timeseries_trunc = my_timeseries.truncate(0, 100000, 2) # doctest: +REMOTE_DATA

Concatenating TimeSeries
========================

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

    >>> eve_ts = ts.TimeSeries(sunpy.data.sample.EVE_TIMESERIES, source='eve') # doctest: +REMOTE_DATA
    >>> goes_ts = ts.TimeSeries(sunpy.data.sample.GOES_XRS_TIMESERIES) # doctest: +REMOTE_DATA
    >>> concatenated_timeseries = goes_ts.concatenate(eve_ts) # doctest: +REMOTE_DATA

Note that the more complex `~sunpy.timeseries.TimeSeriesMetaData` object now has 2 entries and shows details on both:

.. code-block:: python

    >>> concatenated_timeseries.meta # doctest: +REMOTE_DATA
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
