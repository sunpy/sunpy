.. _lightcurve_code_ref:

SunPy lightcurve
================

Overview
--------
One of core classes in SunPy is a LightCurve or timeseries. A number of instruments
are supported through subclasses of the base `~sunpy.lightcurve.LightCurve` class.
To see :ref:`lc-sources` for a list of all of them.

Creating a LightCurve
---------------------
LightCurves can either be creating manually or automatically by downloading
their own data (the most common case). Too create a custom `~sunpy.lightcurve.LightCurve`
see the example in the class documentation below. Subclasses of `~sunpy.lightcurve.LightCurve`
for specific instrument provide their own methods for downloading their data.
For more information see :ref:`lc-sources`.

.. _lc-sources:

Instrument LightCurve Classes
-----------------------------

The generic method to create an instrument-specific LightCurve find the instrument
subclass of interest and follow the following example::

    >>> from sunpy.lightcurve import GOESLightCurve
    >>> from sunpy.time import TimeRange
    >>> tr = TimeRange('2013/07/21', '2013/07/22')
    >>> goes = GOESLightCurve.create(tr)

The `~sunpy.lightcurve.LightCurve.create` method will go off and download the
data needed to populate the instance. The following instrument classes are supported.

.. automodapi:: sunpy.lightcurve
