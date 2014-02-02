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

The lightcurve object currently supports the following data sources

- SDO/EVE
- GOES XRS
- PROBA2/LYRA

For SDO/EVE, please note that only the level OCS and average CSV
files is `currently implemented <http://lasp.colorado.edu/home/eve/data/>`_.

1. Creating a Lightcurve
------------------------

The LightCurve object contains two basic attributes, 'data' and
'meta'.  The 'data' attribute is either a pandas Series object or
a pandas DataFrame (a generalization of the Series object).  These
data objects have very powerful methods for handling data based on
times.  The 'meta' attribute contains any available metadata
information in the original data file.  The type and amount of
metadata available varies from data source to data source since there
is no cross-instrument standard for describing time-series metadata
(unlike FITS files).


2. Inspecting LightCurves
-----------------------
A map contains a number of data-associated attributes. To get a quick look at your map simply
type::

    my_map = sunpy.Map(sunpy.AIA_171_IMAGE)
    my_map
    
This will show a representation of the data as well as some of its associated
attributes. A number of other attributes are also available, for example the date, 
exposure time, map center, xrange, yrange
other::

    map_date = my_map.date
    map_exptime = my_map.exposure_time
    map_center = my_map.center
    map_xrange = my_map.xrange
    map_yrange = my_map.yrange
    
To get a list of all of the attributes check the documentation by typing::

    help(my_map)
    
The meta data for the map is accessed by ::

    header = my_map.meta
    
This references the meta data dictionary with the header information as read from the source
file. 

3. Getting at the data
----------------------
The data in a SunPy LightCurve object is accessible through the data attribute. 



4. Creating a plot of your lightcurve
------------------------------
The SunPy lightcurve object has its own built-in plot methods so that it is easy to
quickly view your lightcurve. To create a plot just type::

    my_lightcurve.peek()
    
This will open a matplotlib plot on your screen.  In addition, to
enable users to modify the plot it is possible to grab the matplotlib
figure object by using the plot() command instead of the show()
command. This makes it possible to use the SunPy lightcurve plot as
the foundation for a more complicated figure.


5. Creating your own lightcurve
------------------------------

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
1440 minutes beginning at the current local time.


