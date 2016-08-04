===========
Lightcurves
===========

Time series data are a fundamental part of many data analysis projects as much
in heliophysics as other areas. SunPy therefore provides a timeseries object to
handle this type of data. This directly superceeds the now depreciated lightcurve
datatype.
Once you've read through this guide check out the :doc:`/code_ref/timeseries`
for a more thorough look at SunPy TimeSeries and to see what data sources it
currently supports.

1. Creating a TimeSeries from a data source
-------------------------------------------

To make things easy, SunPy can download several example files which are used
throughout the docs. These files have names like
`~sunpy.data.sample.EVE_LIGHTCURVE` and `~sunpy.data.sample.GOES_LIGHTCURVE`.
To create the sample `sunpy.timeseries.sources.goes.GOESLightCurve` type the
 following into your interactive Python shell::

    >>> import sunpy
    >>> import sunpy.timeseries
    >>> import sunpy.data.sample
    >>> my_timeseries = sunpy.timeseries.TimeSeries(sunpy.data.sample.GOES_LIGHTCURVE, source='GOES')

If you have not downloaded the data already you should get an error and some
instruction on how to download the sample data.

The variable my_timeseries is a :ref:`timeseries` object. To create one from a
local GOES FITS file try the following::

    >>> my_timeseries = sunpy.timeseries.TimeSeries('/mydirectory/myts.fits', source='GOES')   # doctest: +SKIP

SunPy can automatically detects the source for most FITS files, however timeseries
(and lightcurve) data are stored in a variety of file formats (FITS, txt, csv)
and so it's not always possible to detect the source. For this reason it's good
practice to explicitly select the source for the file.

1.1 Creating a TimeSeries from multiple files
---------------------------------------------

You can create a single time series from multiple files for a give source using
the keyword argument concatenate=True, such as:

    >>> my_timeseries = sunpy.timeseries.TimeSeries('/mydirectory/myts1.fits', '/mydirectory/myts2.fits', source='GOES', concatenate=True)   # doctest: +SKIP

Note these must all be from the same source/instrument. The concatenation method
can be used to make a time series from multiple time time series from different
sources.

2. Creating Custom Time Series
-----------------------------

It's also easy to create a timeseries using custom data in a `~pandas.DataFrame`.

    >>> from sunpy.timeseries import TimeSeries
    >>> my_timeseries = sunpy.timeseries.TimeSeries(DataFrame, {})

The first line imports the lightcurve object.
On the second line empty dictionary is used for metadata, the units dictionary
is also omitted, so the TimeSeries isn't unit aware in this case.
A more comprehensive example is:

    >>> from sunpy.timeseries import TimeSeries
    >>> # To generate some data and the corrisponding dates
    >>> base = datetime.datetime.today()
    >>> dates = [base - datetime.timedelta(minutes=x) for x in range(0, 24 * 60)]
    >>> intensity = np.sin(np.arange(0, 12 * np.pi, ((12 * np.pi) / (24*60))))
    >>> # Create the data DataFrame, header MetaDict and units OrderedDict
    >>> data = DataFrame(intensity, index=dates, columns=['intensity'])
    >>> units = OrderedDict([('intensity', u.W/u.m**2)])
    >>> meta = MetaDict({'key':'value'})
    >>> # Create the time series
    >>> ts_custom = sunpy.timeseries.TimeSeries(data, meta, units)

TimeSeries can also be created from an `~astropy.table`:
#################

Or using a `~numpy.array`:
#################



A LightCurve object must be supplied with some data when it is
created.  The data can either be in your current Python session, in a
local file, or in a remote file.  Let's create some fake data and pass
it into a LightCurve object: ::

    >>> from sunpy.lightcurve import LightCurve
    >>> light_curve = LightCurve.create({"param1": range(24 * 60)})

The first line imports the lightcurve object.  Let's look at the
argument in LightCurve.create.  The argument is a dictionary that
contains a single entry with key "param1" with a value of a list of
1440 entries (from 0 to 1439) - these are our 'fake data'
measurements.  Since no other times are provided, a default set of
times are provided.  You can provide your own times very simply using
the 'index' keyword, as is shown below: ::

    >>> import datetime
    >>> base = datetime.datetime.today()
    >>> dates = [base - datetime.timedelta(minutes=x) for x in range(0, 24 * 60)]
    >>> light_curve = LightCurve.create({"param1": range(24 * 60)}, index=dates)

This gives the measurements "param1" a set of times, in this case,
1440 minutes beginning at the current local time.  Under the hood,
this has created a pandas `~pandas.DataFrame` object with a column name "param1",
with an index of times.





1. Creating a TimeSeries from a data source
-------------------------------------------

To create a `~sunpy.timeseries.TimeSeries` object from one of the supported data
sources you will need to use the
import the object into your session.  Unlike the `~sunpy.map.GenericMap`
object and its instrument-specific subclasses, the instrument sub-classes
of `~sunpy.lightcurve.LightCurve` provide a way to download their own data on
creating. The following example creates a `~sunpy.lightcurve.GOESLightCurve`
for the specified time range: ::

    >>> from sunpy.lightcurve import GOESLightCurve
    >>> from sunpy.time import TimeRange
    >>> tr = TimeRange('2013/07/21', '2013/07/22')
    >>> goes = GOESLightCurve.create(tr)

The `~sunpy.lightcurve.GOESLightCurve` will go off and download the data that is
needed and therefore requires an internet connection.

2. Creating Custom Lightcurve
-----------------------------
It is also very easy to create lightcurves using custom data.

    >>> from sunpy.lightcurve import LightCurve
    >>> light_curve = LightCurve.create({"param1": range(24*60)})

Within `~sunpy.lightcurve.LightCurve.create`, we have a dictionary that contains a single entry with key
``param1`` containing a list of 1440 entries (0-1439). As there are no times provided,
so a default set of times are generated.

A LightCurve object must be supplied with some data when it is
created.  The data can either be in your current Python session, in a
local file, or in a remote file.  Let's create some fake data and pass
it into a LightCurve object: ::

    >>> from sunpy.lightcurve import LightCurve
    >>> light_curve = LightCurve.create({"param1": range(24 * 60)})

The first line imports the lightcurve object.  Let's look at the
argument in LightCurve.create.  The argument is a dictionary that
contains a single entry with key "param1" with a value of a list of
1440 entries (from 0 to 1439) - these are our 'fake data'
measurements.  Since no other times are provided, a default set of
times are provided.  You can provide your own times very simply using
the 'index' keyword, as is shown below: ::

    >>> import datetime
    >>> base = datetime.datetime.today()
    >>> dates = [base - datetime.timedelta(minutes=x) for x in range(0, 24 * 60)]
    >>> light_curve = LightCurve.create({"param1": range(24 * 60)}, index=dates)

This gives the measurements "param1" a set of times, in this case,
1440 minutes beginning at the current local time.  Under the hood,
this has created a pandas `~pandas.DataFrame` object with a column name "param1",
with an index of times.

3. Inspecting maps & Getting at the data
----------------------------------------
A lightcurve holds both data as well as meta data. The meta data for the lightcurve
 is accessed by ::

    >>> header = goes.meta

This references the meta data dictionary with the header information as read
from the source file. A word of caution, many data sources provide little to no
meta data so this variable might be empty.

The data in a SunPy `~sunpy.lightcurve.LightCurve` object is accessible through the
`~sunpy.map.LightCurve.data` attribute.  The data is implemented as a
Pandas `~pandas.DataFrame`, so to get a look at what data you have available ::

    >>> goes.data

You can also get a quick overview of what data you have available like so: ::

    >>> goes.data.info()

LightCurves are columnar data so to get at a particular datum you need to
first index the column then the element you want. To get the names of the
available columns: ::

    >>> goes.data.columns

So you can access the 0th element in the column `xrsa` with: ::

    >>> goes.data['xrsa'][0]

You can also grab all of the data at a particular time: ::

    >>> goes.data['xrsa']['2013-07-21 23:59']

This will return a list of entries with times that match the accuracy of the time
you provide. Finally if you want to get at the x or y values: ::

    >>> x = goes.data.index
    >>> y = goes.data.values

You can read more about indexing at the `pandas documentation website
<http://pandas.pydata.org/pandas-docs/stable/>`_.

3. Plotting
-------------------------------------

The SunPy LightCurve object has its own built-in plot methods so that
it is easy to quickly view your lightcurve. To create a plot just
type:

.. plot::
    :include-source:

    from sunpy.lightcurve import GOESLightCurve
    from sunpy.data.sample import GOES_LIGHTCURVE
    goes = GOESLightCurve.create(GOES_LIGHTCURVE)
    goes.peek()

This will open a matplotlib plot on your screen. The `~sunpy.lightcurve.LightCurve.peek()`
function provides a custom view on the data while `~sunpy.lightcurve.LightCurve.plot()`
provides a more generic plot.

In addition, to enable users to modify the plot it is possible to grab the
matplotlib axes object by using the `~sunpy.lightcurve.LightCurve.plot()` command.
This makes it possible to use the SunPy plot as the foundation for a
more complicated figure. For a bit more information about this and some
examples see :ref:`plotting`. Here is one a more complicated example
which makes use of this methodology.

.. plot::
    :include-source:

    from sunpy.lightcurve import GOESLightCurve
    from sunpy.data.sample import GOES_LIGHTCURVE
    goes = GOESLightCurve.create(GOES_LIGHTCURVE)
    fig = plt.figure()
    ax = goes.plot()
    ax.set_ylim(1e-10, 25e-8)
    ax.set_title('My Plot')
    ax.set_ylabel('Watts')
    plt.show()

Here is another more advanced example. Click the source link to see the that
generated this plot.

.. plot::

    from sunpy.lightcurve import GOESLightCurve
    from sunpy.data.sample import GOES_LIGHTCURVE
    import matplotlib.pyplot as plt
    from datetime import datetime
    from astropy.time import Time

    year = 2012
    month = 6
    day = 1
    time_ranges = ((datetime(year,month,day,17,0,0), datetime(year,month,day,18,0,0)),
                   (datetime(year,month,day,19,5,0), datetime(year,month,day,19,50,0)),
                   (datetime(year,month,day,20,20,0), datetime(year,month,day,21,10,0)))

    goes = GOESLightCurve.create(GOES_LIGHTCURVE)

    plt.figure()
    plt.subplot(211)
    ax = goes.data['xrsb'].plot(color='r')
    plt.ylabel(r'1 - 8 $\AA$')
    plt.title(goes.meta.get('TELESCOP'))
    for time_range in time_ranges:
        plt.axvspan(time_range[0], time_range[1], facecolor='gray', alpha=0.5, label='fit')
    ax.set_xticklabels([])
    plt.subplot(212)
    goes.data['xrsa'].plot()
    plt.ylabel(r'0.5 - 4 $\AA$')
    plt.xlabel(goes.data.index[0].to_datetime().strftime('%Y-%m-%d %H:%m') + ' [UT]')
    plt.show()
