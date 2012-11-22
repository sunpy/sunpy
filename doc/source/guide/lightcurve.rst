===========
Lightcurves
===========

Time series data are a fundamental part of many data analysis projects as much in 
heliophysics as other areas. SunPy therefore provides a lightcurve object to deal 
with this data type. Thankfully most of the heavy lifting in this area has already been
done for us. The lightcurve object makes use of the `pandas <http://pandas.pydata.org/`_
python module. Pandas is a high-quality and optimized module which is in use in a wide 
variety of academic and commercial fields, including Finance, Neuroscience, Economics, 
Statistics, Advertising, and Web Analytics. The lightcurve object is essentially a wrapper
around a pandas dataframe object which holds some of the meta-data from the original 
data source. In this tutorial we provide a quick introduction to 
the lightcurve object and pandas. We highly recommend any user of the lightcurve object 
take a look at the great `pandas documentation <http://pandas.pydata.org/pandas-docs/stable/`_.
for more information.

============
Data Support
============

The lightcurve object currently supports the following data sources

- SDO/EVE
- GOES XRS
- PROBA2/LYRA

1. Creating a Lightcurve
------------------------

A LightCurve object consists of two parts - times, and measurements which were taken at
those times.

A LightCurve object must be supplied with some data when it is created.  The data
can either be in your current Python session, in a local file, or in a remote file.
Let's create some fake data and pass it into a LightCurve object.

    >>> from sunpy.lightcurve import LightCurve
    >>> import datetime
    >>> base = datetime.datetime.today()
    >>> dates = [base - datetime.timedelta(minutes=x) for x in range(0, 24 * 60)]
    >>> light_curve = LightCurve.create({"param1": range(24 * 60)}, index=dates)
    >>> light_curve.show()

The first two lines import sunpy and the datetime module.  The next two lines create a list
of dates.  Let's look at the argument and keyword in LightCurve.create.  The argument is
a dictionary that contains a single entry with key "param1" with a value of a list of
1440 entries (from 0 to 1439).  The keyword "index" passes the dates created earlier.  The
