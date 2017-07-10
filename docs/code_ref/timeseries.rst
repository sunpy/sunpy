.. _timeseries_code_ref:

SunPy timeseries
================

Overview
--------
One of the core classes in SunPy is a timeseries. A number of instruments
are supported through subclasses of the base `~sunpy.timeseries.GenericTimeSeries` class.
To see :ref:`ts-sources` for a list of all of them.

Creating a TimeSeries
---------------------
TimeSeries can either be created manually or from source files using the factory.
To create a custom `~sunpy.timeseries.GenericTimeSeries`
see the example in the class documentation below. Subclasses of `~sunpy.timeseries.GenericTimeSeries`
for specific instruments provide their own methods for opening files for their data.
For more information see :ref:`ts-sources`.

.. _ts-sources:

Instrument TimeSeries Classes
-----------------------------

The generic method to create an instrument-specific TimeSeries is to call the `~sunpy.timeseries.TimeSeries`
factory on a file from that dataset and pass the source keyword argument. Some sources use
FITS files and these may define the source class themselves, in this case the TimeSeries
factory can determine the source implicitly, but it's good practice to explicitly state it.
The following example shows the factory loading a sample file::

    >>> import sunpy.timeseries as ts
    >>> import sunpy.data.sample
    >>> goes = ts.TimeSeries(sunpy.data.sample.GOES_XRS_TIMESERIES, source='XRS')

The `~sunpy.timeseries.TimeSeries` factory will load the file and create the timeseries
instance. The following instrument classes are supported.

.. automodapi:: sunpy.timeseries
