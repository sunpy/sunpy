===========
Lightcurves
===========

Time series data are a fundamental part of many data analysis projects as much in 
heliophysics as other areas. SunPy therefore provides a lightcurve object to deal 
with this data type. Thankfully most of the heavy lifting in this area has already been
done for us. The lightcurve object makes use of the `pandas <http://pandas.pydata.org/>`_
python module. Pandas is a high-quality and optimized module which is in use in a wide 
variety of academic and commercial fields, including Finance, Neuroscience, Economics, 
Statistics, Advertising, and Web Analytics. The lightcurve object is essentially a wrapper
around a pandas dataframe object which holds some of the meta-data from the original 
data source. In this tutorial we provide a quick introduction to 
the lightcurve object and pandas. We highly recommend any user of the lightcurve object 
take a look at the great `pandas documentation <http://pandas.pydata.org/pandas-docs/stable/>`_.
for more information.

Data Support
------------

The lightcurve object currently supports the following data sources:

- SDO/EVE
- GOES XRS
- PROBA2/LYRA

For SDO/EVE, please note that only the level OCS and average CSV
files is `currently implemented <http://lasp.colorado.edu/home/eve/data/>`_.

1. Creating a Lightcurve from a data source
-------------------------------------------

To create a LightCurve object from one of the supported data sources,
import the object into your Python session.  The example below gets
LYRA data for date indicated:

    >>> from sunpy.lightcurve import LYRALightCurve
    >>> ly = LYRALightCurve.create('2013/07/21')

The LYRA LightCurve object parses the input date and downloads the
relevant data from the PROBA2/LYRA webpage.

LightCurve objects for supported data sources can also be created by
passing in a local file path and file name or by specifying the URL
where the data is located.


2.  Getting at the data
-----------------------

The data in a SunPy LightCurve object is accessible through the data
attribute.  The header is accessible through the meta attribute.  Note
that, unlike FITS files, there is no cross-instrument standard
governing what should be kept in a time-series arising from a
particular instrument.

The data attribute stores the time-series measurements as a pandas
DataFrame object.  For example,

    >>> ly.data
    >>> <class 'pandas.core.frame.DataFrame'>
    ... DatetimeIndex: 1723918 entries, 2013-07-21 00:00:00.129000 to 2013-07-22 00:00:00.033000
    ... Columns: 4 entries, CHANNEL1 to CHANNEL4
    ... dtypes: float64(4)

Pandas has a lot of extremely useful methods for manipulating
time-series information.  Please consult the `pandas documentation
<http://pandas.pydata.org/pandas-docs/stable/>`_ to find out how to
access the data inside a pandas DataFrame object.


3. Creating a plot of your lightcurve
-------------------------------------

The SunPy LightCurve object has its own built-in plot methods so that
it is easy to quickly view your lightcurve. To create a plot just
type:

    >>> my_lightcurve.peek()
    
This will open a matplotlib plot on your screen.  In addition, to
enable users to modify the plot, it is possible to grab the matplotlib
figure object by using the plot() command instead of the show()
command. This makes it possible to use the SunPy lightcurve plot as
the foundation for a more complicated figure.


4. Creating your own lightcurve
-------------------------------

A LightCurve object must be supplied with some data when it is
created.  The data can either be in your current Python session, in a
local file, or in a remote file.  Let's create some fake data and pass
it into a LightCurve object.

    >>> from sunpy.lightcurve import LightCurve
    >>> light_curve = LightCurve.create({"param1": range(24 * 60)})

The first line imports the lightcurve object.  Let's look at the
argument in LightCurve.create.  The argument is a dictionary that
contains a single entry with key "param1" with a value of a list of
1440 entries (from 0 to 1439) - these are our 'fake data'
measurements.  Since no other times are provided, a default set of
times are provided.  You can provide your own times very simply using
the 'index' keyword, as is shown below.

    >>> import datetime
    >>> base = datetime.datetime.today()
    >>> dates = [base - datetime.timedelta(minutes=x) for x in range(0, 24 * 60)]
    >>> light_curve = LightCurve.create({"param1": range(24 * 60)}, index=dates)

This gives the measurements "param1" a set of times, in this case,
1440 minutes beginning at the current local time.  Under the hood,
this has created a pandas DataFrame object with a colum name "param1",
with an index of times.


