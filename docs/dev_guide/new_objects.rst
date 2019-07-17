.. _new_maps_ts_etc:

************************************************
Creating new SunPy Subclasses (Maps, TimeSeries)
************************************************

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
