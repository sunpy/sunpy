***********
TimeSeries
***********

Time series data are a fundamental part of many data analysis projects
in heliophysics as well as other areas. SunPy therefore provides a TimeSeries
object to handle this type of data. This new object supersedes the now
deprecated Lightcurve datatype. Once you've read through this guide check out
the :doc:`/code_ref/timeseries` for a more thorough look at SunPy TimeSeries
and to see what data sources it currently supports.

Sections 1 and 2 below describe how to create TimeSeries objects.  Section 3
describes how to inspect the TimeSeries object to examine the data itself and
the metadata.  Section 4 describes how to plot a TimeSeries object, Section 5
describes how a TimeSeries object can be manipulated, and Section 6 describes
in detail the TimeSeries metadata.

.. warning::

   The TimeSeries supersedes the previous LightCurve object but does not
   implement data download methods. To download TimeSeries data files use
   `sunpy.net`.

1. Creating a TimeSeries from an observational data source
==========================================================

A TimeSeries object can be created from local files.  For convenience, SunPy can
download several example timeseries of observational data. These files have names like
``sunpy.data.sample.EVE_TIMESERIES`` and ``sunpy.data.sample.GOES_XRS_TIMESERIES``.
To create the sample `sunpy.timeseries.sources.goes.XRSTimeSeries` type the
following into your interactive Python shell: ::

    >>> import sunpy.timeseries as ts
    >>> import sunpy.data.sample  # doctest: +REMOTE_DATA
    >>> my_timeseries = ts.TimeSeries(sunpy.data.sample.GOES_XRS_TIMESERIES, source='XRS')  # doctest: +REMOTE_DATA

.. doctest-skip-all

This is calling the `~sunpy.timeseries.TimeSeries` factory to create a time
series from a GOES XRS FITS file. The TimeSeries factory uses `sunpy.io.fits` to
read the FITS file. Note that if you have not downloaded the data already you
should get an error and some instruction on how to download the sample data.

The variable ``my_timeseries`` is a :ref:`timeseries` object. To create one from
a local GOES/XRS FITS file try the following: ::

    >>> my_timeseries = ts.TimeSeries('/mydirectory/myts.fits', source='XRS')   # doctest: +SKIP

SunPy will attempt to detect automatically the instrument source for most FITS
files. However timeseries data are stored in a variety of file types (FITS, txt,
csv) and formats, and so it is not always possible to detect the source. For
this reason, it is good practice to explicitly state the source for the file.
SunPy ships with a number of known instrumental sources.  If you would like
SunPy to include another instrumental source please follow the `SunPy
contribution guide <https://sunpy.org/contribute/>`.

The `~sunpy.timeseries.TimeSeries` factory has the ability to make a list of
TimeSeries objects using a list of filepaths, a folder or a glob, for example: ::

    >>> my_ts_list = ts.TimeSeries('filepath1', 'filepath2', source='XRS')   # doctest: +SKIP
    >>> my_ts_list = ts.TimeSeries('/goesdirectory/', source='XRS')   # doctest: +SKIP
    >>> my_ts_list = ts.TimeSeries(glob, source='XRS')   # doctest: +SKIP

Note that this functionality will only work with files from the same single
source, generating a source specific child of the `~sunpy.timeseries.GenericTimeSeries`
class such as the `~sunpy.timeseries.sources.goes.XRSTimeSeries` above. For this
reason, all the files should be from that same source for the `~sunpy.timeseries.TimeSeries`
factory to work as intended.

1.1 Creating a Single TimeSeries from Multiple Files
====================================================

You can create a single time series from multiple files for a given source using
the keyword argument ``concatenate=True``, such as:

    >>> my_timeseries = ts.TimeSeries('/mydirectory/myts1.fits', '/mydirectory/myts2.fits', source='XRS', concatenate=True)  # doctest: +SKIP

Note these must all be from the same source if using
`~sunpy.timeseries.GenericTimeSeries.concatenate` from within the TimeSeries
factory. However, if time series `~sunpy.timeseries.GenericTimeSeries.concatenate` method
can be used to make a single time series from multiple TimeSeries from different
sources once they are already in the form of TimeSeries objects.

2. Creating Custom TimeSeries
=============================

Sometimes you will have data that you want to create into a TimeSeries. You can
use the factory to create a `~sunpy.timeseries.GenericTimeSeries`
from a variety of data sources currently including `pandas.DataFrame` and
`astropy.table.Table`.

2.1 Creating a TimeSeries from a Pandas DataFrame
=================================================

A TimeSeries object must be supplied with some data when it is
created.  The data can either be in your current Python session, in a
local file, or in a remote file.  Let's create some data and pass
it into a TimeSeries object: ::

    >>> import numpy as np
    >>> intensity = np.sin(np.arange(0, 12 * np.pi, ((12 * np.pi) / (24*60))))

The first line imports the numpy module used to create and store the data.
The second line creates a basic numpy array of values representing a sine wave.
We can use this array along with a suitable time storing object (such as Astropy
`~astropy.time` or a list of `datetime` objects) to make a Pandas
`~pandas.DataFrame`.  A suitable list of times must contain the same
number of values as the data, this can be created using: ::

    >>> import datetime
    >>> base = datetime.datetime.today()
    >>> times = [base - datetime.timedelta(minutes=x) for x in range(24*60, 0, -1)]

The Pandas `~pandas.DataFrame` will use the dates list as the index: ::

    >>> from pandas import DataFrame
    >>> data = DataFrame(intensity, index=times, columns=['intensity'])

This `~pandas.DataFrame` can then be used to construct a TimeSeries: ::

    >>> import sunpy.timeseries as ts
    >>> ts_custom = ts.TimeSeries(data)

Furthermore we could specify the metadata/header and units of this time series
by sending them as arguments to the factory: ::

    >>> from collections import OrderedDict
    >>> import astropy.units as u

    >>> meta = OrderedDict({'key':'value'})
    >>> units = OrderedDict([('intensity', u.W/u.m**2)])
    >>> ts_custom = ts.TimeSeries(data, meta, units)

2.2 Creating Custom TimeSeries from an Astropy Table
====================================================

A Pandas `~pandas.DataFrame` is the underlying object used to store
the data within a TimeSeries, so the above example is the most lightweight to
create a custom TimeSeries, but being scientific data it will often be more
convenient to use an Astropy `~astropy.table.Table` and let the factory
convert this.  An advantage of this method is it allows you to include metadata
and Astropy `~astropy.units.quantity.Quantity` values, which are both supported
in tables, without additional arguments.  For example: ::

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

Note that due to the properties of the `~astropy.time.Time` object, this will be
a mixin column which since it is a single object, limits the versatility of
the `~astropy.table.Table` a little. For more on mixin columns see the `Astropy
docs <https://docs.astropy.org/en/stable/table/mixin_columns.html>`_.  The units
will be taken from the table quantities for each column, the metadata will
simply be the table.meta dictionary.  You can also explicitly add metadata and
units, these will be added to the relevant dictionaries using the dictionary
update method, with the explicit user-given values taking precedence.

    >>> from sunpy.util.metadata import MetaDict
    >>> from collections import OrderedDict
    >>> import astropy.units as u

    >>> meta = MetaDict({'key':'value'})
    >>> units = OrderedDict([('intensity', u.W/u.m**2)])
    >>> ts_table = ts.TimeSeries(table, meta, units)


3. Inspecting TimeSeries & Getting at the Data
===============================================

A time series holds both data as well as meta data and units data. The meta data
for the time series is accessed by: ::

    >>> header = my_timeseries.meta

This references the `~sunpy.timeseries.TimeSeriesMetaData` object with
the header information as read from the source files. A word of caution: many
data sources provide little to no meta data so this variable might be empty.
The meta data is described in more detail later in this guide. Similarly there
are properties for getting `~sunpy.timeseries.GenericTimeSeries.columns`
as a list of strings, `~sunpy.timeseries.GenericTimeSeries.get_index`
values and `~sunpy.timeseries.GenericTimeSeries.time_range` of
the data.  The actual data in a SunPy TimeSeries object is accessible through
the `~sunpy.timeseries.GenericTimeSeries.data` attribute.  The
data is implemented as a Pandas `~pandas.DataFrame`, so to get a look at what
data you have available use: ::

    >>> my_timeseries.data  # doctest: +SKIP

You can also get a quick overview of that data using: ::

    >>> my_timeseries.data.info()
    <class 'pandas.core.frame.DataFrame'>
    DatetimeIndex: 42177 entries, 2011-06-06 23:59:59.961999 to 2011-06-07 23:59:57.631999
    Data columns (total 2 columns):
    xrsa    42177 non-null float32
    xrsb    42177 non-null float32
    dtypes: float32(2)
    memory usage: 659.0 KB

Time series are columnar data so to get at a particular datum you need to
first index the column, then the element you want. To get the names of the
available columns: ::

    >>> my_timeseries.data.columns
    Index(['xrsa', 'xrsb'], dtype='object')

You can access the 0th element in the column ``xrsa`` with: ::

    >>> my_timeseries.data['xrsa'][0]
    1e-09

You can also grab all of the data at a particular time: ::

    >>> my_timeseries.data['xrsa']['2011-06-07 00:00:02.008999']
    1e-09

This will return a list of entries with times that match the accuracy of the time
you provide. You can consider the data as x or y values: ::

    >>> x = my_timeseries.data.index
    >>> y = my_timeseries.data.values

You can read more about indexing at the `pandas documentation website
<https://pandas.pydata.org/pandas-docs/stable/>`_.

A TimeSeries can also return an Astropy `~astropy.units.quantity.Quantity` for a
given column using the `~sunpy.timeseries.GenericTimeSeries.quantity`
method, this uses the values stored in the data and units stored in the units
dictionary to determine the `~astropy.units.quantity.Quantity`: ::

    >>> quantity = my_timeseries.quantity('xrsa')

4. Plotting
===========

The SunPy TimeSeries object has its own built-in plot methods so that
it is easy to quickly view your time series. To create a plot just
type:

.. plot::
    :include-source:

    import sunpy.timeseries as ts
    import sunpy.data.sample
    ts_plot = ts.TimeSeries(sunpy.data.sample.GOES_XRS_TIMESERIES, source='XRS')
    fig = ts_plot.peek()

This will open a matplotlib plot on your screen. The `~sunpy.timeseries.GenericTimeSeries.peek`
method provides a view on data customised for each source while `~sunpy.timeseries.GenericTimeSeries.plot`
provides a more general plot.  Note that `~sunpy.timeseries.GenericTimeSeries.peek`
returns a `matplotlib.figure.Figure` object, if you want to save this to a PNG
file you can use the ``savefig`` method:

    >>> fig.savefig('figure.png')  # doctest: +SKIP

In addition, to enable users to modify the plot it is possible to grab the
matplotlib axes object by using the `~sunpy.timeseries.GenericTimeSeries.plot`
command.  This makes it possible to use the SunPy plot as the foundation for a
more complicated figure. For a more information about this and some examples see
:ref:`plotting`.


5 Manipulating TimeSeries
=========================

5.1 Modifying the Data
======================

Since the timeseries data is stored as a Pandas `~pandas.DataFrame`
you can easily modify the data directly using all of the usual Pandas methods:
for example, you can modify a single cells value using: ::

    >>> my_timeseries.data['xrsa'][0] = 0.1

Or similarly using a datetime values (as string or datetime object): ::

    >>> my_timeseries.data['xrsa']['2012-06-01 23:59:45.061999'] = 1

You can even change all the values for a given time: ::

    >>> my_timeseries.data['xrsa']['2012-06-01 00:00'] = 1

Note, you will need to be careful to consider units when modifying the
TimeSeries data directly. For further details about editing Pandas DataFames you
can read the `pandas documentation website <https://pandas.pydata.org/pandas-docs/stable/>`_.

Additionally the TimeSeries provides the `~sunpy.timeseries.GenericTimeSeries.add_column`
method which will either add a new column or update a current column if the
colname is already present. This can take numpy array or preferably an Astropy
`~astropy.units.quantity.Quantity` value.  For example: ::

    >>> values = u.Quantity(my_timeseries.data['xrsa'].values[:-2], my_timeseries.units['xrsa']) * 20.5
    >>> my_timeseries.add_column('new col', values)
    <sunpy.timeseries.sources.goes.XRSTimeSeries object at ...>

Note that the values will be converted into the column units if an Astropy
`~astropy.units.quantity.Quantity` is given. Caution should be taken when adding
a new column because this column won't have any associated MetaData entry,
similarly if you use an array of values it won't add an entry into the units
`~collections.OrderedDict`.

5.2 Truncating a TimeSeries
===========================

It is often useful to truncate an existing TimeSeries object to retain a
specific time range.  This is easily achieved by using the `~sunpy.timeseries.GenericTimeSeries.truncate`
method. For example, to trim our GOES data into a period of interest use: ::

    >>> from sunpy.time import TimeRange
    >>> tr = TimeRange('2012-06-01 05:00','2012-06-01 06:30')
    >>> my_timeseries_trunc = my_timeseries.truncate(tr)

This takes a number of different arguments, such as the start and end dates (as
datetime or string objects) or a `~sunpy.time.TimeRange` as used above. Note
that the truncated TimeSeries will have a truncated `~sunpy.timeseries.TimeSeriesMetaData`
object, which may include dropping metadata entries for data totally cut out
from the TimeSeries.  If you want to truncate using slice-like values you can,
for example taking every 2nd value from 0 to 10000 can be done using: ::

    >>> my_timeseries_trunc = my_timeseries.truncate(0,100000,2)

Caution should be used when removing values from the data manually, the
TimeSeries can't guarantee Astropy units are correctly preserved when you
interact with the data directly.

5.3 Down and Up Sampling a TimeSeries Using Pandas
==================================================

Because the data is stored in a Pandas `~pandas.DataFrame` object you
can manipulate it using normal Pandas methods, such as the `~pandas.DataFrame.resample`
method.  To downsample you can use: ::

    >>> downsampled_dataframe = my_timeseries_trunc.data.resample('10T').mean()

Note, here ``10T`` means sample every 10 minutes and 'mean' is the method used
to combine the data. Alternatively the sum method is often used.
You can also upsample, such as: ::

    >>> upsampled_data = my_timeseries_trunc.data.resample('30S').ffill()

Note, here we upsample to 30 second intervals using ``30S`` and use the pandas
fill-forward method. Alternatively the back-fill method could be used.  Caution
should be used when resampling the data, the TimeSeries can't guarantee Astropy
Units are correctly preserved when you interact with the data directly.

5.4 Concatenating TimeSeries
============================

It's common to want to combine a number of TimeSeries together into a single
TimeSeries.  In the simplest scenario this is to combine data from a single
source over several time ranges, for example if you wanted to combine the daily
GOES data to get a week or more of constant data in one TimeSeries.  This can be
performed using the TimeSeries factory with the ``concatenate=True``
keyword argument: ::

    >>> concatenated_timeseries = sunpy.timeseries.TimeSeries(filepath1, filepath2, source='XRS', concatenate=True)  # doctest: +SKIP

Note, you can list any number of files, or a folder or use a glob to select the
input files to be concatenated.  It is possible to concatenate two TimeSeries
after creating them with the factory using the `~sunpy.timeseries.GenericTimeSeries.concatenate`
method.  For example: ::

    >>> concatenated_timeseries = goes_timeseries_1.concatenate(goes_timeseries_2)  # doctest: +SKIP

This will result in a TimeSeries identical to if you used the factory to create
it in one step.  A limitation of the TimeSeries class is that often it is not
easy to determine the source observatory/instrument of a file, generally
because the file formats used vary depending on the scientific working groups,
thus some sources need to be explicitly stated (as a keyword argument) and so it
is not possible to concatenate files from multiple sources with the factory.
To do this you can still use the `~sunpy.timeseries.GenericTimeSeries.concatenate`
method, which will create a new TimeSeries with all the rows and columns of the
source and concatenated TimeSeries in one: ::

    >>> concatenated_timeseries = goes_timeseries.concatenate(eve_timeseries)  # doctest: +SKIP

Note that the more complex `~sunpy.timeseries.TimeSeriesMetaData`
object now has 2 entries and shows details on both: ::

    >>> concatenated_timeseries.meta  # doctest: +SKIP

The metadata object is described in more detail in the next section.


5.5 Creating an Astropy Table from a TimeSeries
===============================================

If you want to take the data from your TimeSeries and use it as a `~astropy.table.Table`
this can be done using the `~sunpy.timeseries.GenericTimeSeries.to_table`
method.  For example: ::

    >>> table = my_timeseries_trunc.to_table()

Note that this `~astropy.table.Table` will contain a mixin column for
containing the Astropy `~astropy.time.Time` object representing the index,
it will also add the relevant units to the columns. One of the most useful
reasons for doing this is that Astropy `~sunpy.timeseries.GenericTimeSeries.to_table`
objects have some very nice options for viewing the data, including the basic
console view: ::

    >>> table
    <Table length=21089>
                 date               xrsa     xrsb
                                   W / m2   W / m2
            datetime64[ns]        float32  float32
    ----------------------------- ------- ----------
    2011-06-06T23:59:59.961999000     0.1 1.8871e-07
    2011-06-07T00:00:04.058999000   1e-09 1.8609e-07
    2011-06-07T00:00:08.151999000   1e-09 1.8609e-07
    2011-06-07T00:00:12.248999000   1e-09 1.8609e-07
    2011-06-07T00:00:16.344999000   1e-09 1.8084e-07
    2011-06-07T00:00:20.441999000   1e-09 1.8084e-07
    2011-06-07T00:00:24.534999000   1e-09 1.8084e-07
    2011-06-07T00:00:28.631999000   1e-09 1.8346e-07
    2011-06-07T00:00:32.728999000   1e-09 1.8346e-07
                              ...     ...        ...
    2011-06-07T23:59:20.768999000   1e-09  1.651e-07
    2011-06-07T23:59:24.864999000   1e-09 1.5985e-07
    2011-06-07T23:59:28.961999000   1e-09 1.5985e-07
    2011-06-07T23:59:33.058999000   1e-09 1.6248e-07
    2011-06-07T23:59:37.151999000   1e-09 1.6248e-07
    2011-06-07T23:59:41.248999000   1e-09 1.5985e-07
    2011-06-07T23:59:45.344999000   1e-09 1.5723e-07
    2011-06-07T23:59:49.441999000   1e-09 1.6248e-07
    2011-06-07T23:59:53.538999000   1e-09 1.5985e-07
    2011-06-07T23:59:57.631999000   1e-09 1.5985e-07

and the more sophisticated browser view using the `~astropy.table.Table.show_in_browser`
method: ::

    >>> table.show_in_browser(jsviewer=True)  # doctest: +SKIP

For further details about editing Astropy tables you can read the `astropy
documentation website <https://docs.astropy.org/en/stable/table/>`_.


6. A Detailed Look at the Metadata
==================================

TimeSeries store metadata in a `~sunpy.timeseries.TimeSeriesMetaData`
object, this object is designed to be able to store multiple basic `~sunpy.util.metadata.MetaDict`
(case-insensitive ordered dictionary) objects and able to identify the relevant
metadata for a given cell in the data. This enables a single TimeSeries to be
created by combining/concatenating multiple TimeSeries source files together
into one and to keep a reliable track of all the metadata relevant to each cell,
column or row.  The metadata can be accessed by: ::

    >>> meta = my_timeseries.meta

You can easily get an overview of the metadata, this will show you a basic
representation of the metadata entries that are relevant to this TimeSeries. ::

    >>> meta
    |-------------------------------------------------------------------------------------------------|
    |TimeRange                  | Columns         | Meta                                              |
    |-------------------------------------------------------------------------------------------------|
    |2011-06-06 23:59:59.961999 | xrsa            | simple: True                                      |
    |            to             | xrsb            | bitpix: 8                                         |
    |2011-06-07 23:59:57.631999 |                 | naxis: 0                                          |
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

The data within a `~sunpy.timeseries.TimeSeriesMetaData` object is
stored as a list of tuples, each tuple representing the metadata from a source
file or timeseries. The tuple will contain a `~sunpy.time.TimeRange` telling us
which rows the metadata applies to, a list of column name strings for which the
metadata applies to and finally a `~sunpy.util.metadata.MetaDict` object for
storing the key/value pairs of the metadata itself.  Each time a TimeSeries is
concatenated to the original a new set of rows and/or columns will be added to
the `~pandas.DataFrame` and a new entry will be added into the
metadata.  Note that entries are ordered chronologically based on
`~sunpy.time.timerange.TimeRange.start` and generally it's expected that no two
TimeSeries will overlap on both columns and time range.  For example it is not
good practice for alternate row values in a single column to be relevant to
different metadata entries as this would make it impossible to uniquely identify
the metadata relevant to each cell.

If you want the string that's printed then you can use the
`~sunpy.timeseries.TimeSeriesMetaData.to_string` method.  This has the
advantage of having optional keyword arguments that allows you to set the depth
(number of rows for each entry) and width (total number of characters wide)
to better fit your output.  For example: ::

    >>> meta_str = meta.to_string(depth = 20, width=99)

Similar to the TimeSeries, the metadata has some properties for convenient
access to the global metadata details, including
`~sunpy.timeseries.TimeSeriesMetaData.columns` as a list of
strings, `~sunpy.timeseries.TimeSeriesMetaData.index` values
and `~sunpy.timeseries.TimeSeriesMetaData.time_range` of the data.
Beyond this, there are properties to get lists of details for all the entries in
the `~sunpy.timeseries.TimeSeriesMetaData` object, including
`~sunpy.timeseries.TimeSeriesMetaData.timeranges`,
`~sunpy.timeseries.TimeSeriesMetaData.columns` (as a list of string
column names) and `~sunpy.timeseries.TimeSeriesMetaData.metas`.
Similar to TimeSeries objects you can `~sunpy.timeseries.TimeSeriesMetaData.truncate`
and `~sunpy.timeseries.TimeSeriesMetaData.concatenate` `~sunpy.timeseries.TimeSeriesMetaData`
objects, but generally you won't need to do this as it is done automatically
when actioned on the TimeSeries.
Note that when truncating a `~sunpy.timeseries.TimeSeriesMetaData`
object you will remove any entries outside of the given `~sunpy.time.TimeRange`.
You can also `~sunpy.timeseries.TimeSeriesMetaData.append` a new entry
(as a tuple or list), which will add the entry in the correct chronological
position.  It is frequently necessary to locate the metadata for a given column,
row or cell which can be uniquely identified by both, to do this you can use the
`~sunpy.timeseries.TimeSeriesMetaData.find` method, by adding colname
and/or time/row keyword arguments you get a `~sunpy.timeseries.TimeSeriesMetaData`
object returned which contains only the relevant entries. You can then use the
`~sunpy.timeseries.TimeSeriesMetaData.metas` property to get a list of
just the relevant `~sunpy.util.metadata.MetaDict` objects.  For example: ::

    >>> tsmd_return = my_timeseries.meta.find(colname='xrsa', time='2012-06-01 00:00:33.904999')
    >>> tsmd_return.metas
    []

Note, the colname and time filters are optional, but omitting both filters just
returns an identical `~sunpy.timeseries.TimeSeriesMetaData` object to
the TimeSeries original. A common use case for the metadata is to find out the
instrument/s that gathered the data and in this case you can use the
`~sunpy.timeseries.TimeSeriesMetaData.get` method.  This method takes a
single key string or list of key strings with the optional filters and will
search for any matching values. This method returns another `~sunpy.timeseries.TimeSeriesMetaData`
object, but removes all unwanted key/value pairs.  The result can be converted
into a simple list of strings using the `~sunpy.timeseries.TimeSeriesMetaData.values`
method: ::

    >>> tsmd_return = my_timeseries.meta.get('telescop', colname='xrsa')
    >>> tsmd_return.values()
    ['GOES 15']

Note `~sunpy.timeseries.TimeSeriesMetaData.values` removes duplicate
strings and sorts the returned list.  You can update the values for these
entries efficiently using the `~sunpy.timeseries.TimeSeriesMetaData.update`
method which takes a dictionary argument and updates the values to each of the
dictionaries that match the given colname and time filters, for example: ::

    >>> my_timeseries.meta.update({'telescop': 'G15'}, colname='xrsa', overwrite=True)

Here we have to specify the overwrite=False keyword parameter to allow us to
overwrite values for keys already present in the `~sunpy.util.metadata.MetaDict`
objects, this helps protect the integrity of the original metadata and without
this set (or with it set to False) you can still add new key/value pairs.
Note that the `~sunpy.util.metadata.MetaDict` objects are both case-insensitive
for key strings and have ordered entries, where possible the order is preserved
when updating values.
