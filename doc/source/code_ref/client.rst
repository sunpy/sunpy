.. _dataretriever_code_ref:

SunPy GenericClient
===================

Overview
--------
One of core classes in SunPy is a GenericClient or dataretriver . A number of instruments
are supported through subclasses of the base `~sunpy.net.dataretriever.GenericClient` class.
See :ref:`dataretriever-sources` for a list of all of them.

Creating a Client and querying for data
---------------------------------------
Clients can be manually created and used for downloading
lightcurve instrument data. To create a custom `~sunpy.net.dataretiever.GenericClient`
see the examples in the instrument/sources documentation below. Subclasses of `~sunpy.net.dataretiever.GenericClient`
for specific instrument provide their own methods for downloading their data.
For more information see :ref:`dataretriever-sources`.

.. _dataretriever-sources:

Instrument Client Classes
-----------------------------

The generic method to create an instrument-specific Client is to find the instrument
subclass of interest and follow the following example::

	    >>> from sunpy.time.timerange import TimeRange
	    >>> from sunpy.net.vso.attrs import Time, Instrument
	    >>> import sunpy.net.dataretriever.sources.stereo as stereo
	    >>> LCClient = stereo.PLASTICClient()

	    >>> qr1 = LCClient.query(Time(TimeRange('2012/11/27', '2012/11/27')), Instrument('stereo/plastic'))

Example Above returns a response object which is further used to download the corresponding data file as shown below::

	    >>> res = LCClient.get(qr1)
	    >>> download_list = res.wait()

Doing so will download the data files corresponding to input timerange and instrument and Default values for other client specific parameters. The following instrument classes are supported.

.. automodapi:: sunpy.net.dataretriever

