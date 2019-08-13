.. _new_maps_ts_etc:

************************************************
Creating new SunPy Subclasses (Maps, TimeSeries)
************************************************

Documenting new data sources
----------------------------
Each subclass of `~sunpy.map.GenericMap` or `~sunpy.timeseries.TimeSeries`
must provide a detailed docstring providing an overview of the data source that
the object represents. In order to maintain consistency and completeness the
following information must be provided by a data source docstring, if available,
and preferably in the following order:

* the name of the mission and instrument and the institution that built it
* a short description of the instrument (e.g. Cassegrain reflector, Wolter-1 grazing incidence x-ray, coronagraph) including the type of detector
* a description of the platform (e.g. satellite in 28 deg inclined orbit, a telescope on the summit of Mauna Kea in Hawaii)
* a description of the primary purpose or science goals of the instrument.
* a list of all wavelength(s) or passbands in appropriate units
* a description of the emission processes which dominate in those passbands
* the field of view and resolution (e.g. angular resolution)
* the measurement cadence
* a description of the operational concept (e.g. operates 24/7, observes from 7 am to 5 pm UT)
* the start and end of the data set

In addition, a reference section must be provided with links to the following
resources, if available,

* the mission web page
* the instrument web page
* relevant wikipedia page(s)
* relevant user guide(s)
* reference to the mission paper
* reference to the instrument paper
* information to interpret metadata keywords such as FITS header reference
* link to the data archive

Writing a new Instrument Map Class
==================================

Any subclass of `~sunpy.map.GenericMap` which defines a method named `~sunpy.map.GenericMap.is_datasource_for` will automatically be registered with the `Map <sunpy.map.map_factory.MapFactory>` factory.
The ``is_datasource_for`` method describes the form of the data and metadata for which the `~sunpy.map.GenericMap` subclass is valid.
For example it might check the value of the ``INSTRUMENT`` key in the metadata dictionary.
This makes it straightforward to define your own `~sunpy.map.GenericMap` subclass for a new instrument or a custom data source like simulated data.
These classes only have to be imported for this to work, as demonstrated by the following example.

.. code-block:: python

    import sunpy.map
    class FutureMap(sunpy.map.GenericMap):

        def __init__(self, data, header, **kwargs):

            super(FutureMap, self).__init__(data, header, **kwargs)

            # Any Future Instrument specific keyword manipulation

       # Specify a classmethod that determines if the data-header pair matches
       # the new instrument
       @classmethod
       def is_datasource_for(cls, data, header, **kwargs):
            """Determines if header corresponds to an AIA image"""
            return header.get('instrume', '').startswith('FUTURESCOPE')


This class will now be available through the `Map <sunpy.map.map_factory.MapFactory>` factory as long as this class has been defined, i.e. imported into the current session.

If you do not want to create a method named ``is_datasource_for`` you can manually register your class and matching method using the following method

.. code-block:: python

    import sunpy.map

    sunpy.map.Map.register(FutureMap, FutureMap.some_matching_method)

Writing a new Instrument TimeSeries Class
=========================================

To be written.
