.. _timeseries_code_ref:

Timeseries (`sunpy.timeseries`)
*******************************
One of the core classes in sunpy is a timeseries.
A number of instruments are supported through subclasses of the base `~sunpy.timeseries.GenericTimeSeries` class.
See :ref:`ts-sources` for a list of them.

A timeseries can be created by calling `~sunpy.timeseries.TimeSeries`.

.. _ts-sources:

Instrument TimeSeries Classes
=============================
The generic method to create an instrument-specific TimeSeries is to call `~sunpy.timeseries.TimeSeries` with a file path and the instrument-specific source keyword argument.
In some cases the source can be determined automatically if a FITS file is being loaded.

The following example shows the factory loading a sample file::

    >>> import sunpy.timeseries as ts
    >>> import sunpy.data.sample  # doctest: +REMOTE_DATA
    >>> goes = ts.TimeSeries(sunpy.data.sample.GOES_XRS_TIMESERIES, source='XRS')  # doctest: +REMOTE_DATA

The `~sunpy.timeseries.TimeSeries` factory will load the file and create the timeseries instance.
The following instrument classes are supported:

.. automodapi:: sunpy.timeseries
    :no-inheritance-diagram:
    :include-all-objects:

.. automodapi:: sunpy.timeseries.sources

CDF files
=========
`~sunpy.timeseries.GenericTimeSeries` can load a single CDF file, or a list of CDF files if ``concatenate=True`` is passed.

Units
-----
The physical units of different columns in CDF files do not conform to a standard that `astropy.units` understands.
sunpy internally stores a set of common mappings from unit strings to `~astropy.units.Unit`, but you may see a warning about unrecognised unit strings when reading a CDF file.
To register the correct unit definition :func:`astropy.units.add_enabled_units` can be used.
For example, to register 'deg K' as representing Kelvin and '#/cc' as 1/cm^3::

  >>> import astropy.units as u
  >>> _ = u.add_enabled_units([u.def_unit('deg K', represents=u.K), u.def_unit('#/cc', represents=u.cm**-3)])
