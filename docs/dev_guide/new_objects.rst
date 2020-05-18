.. _new_maps_ts_etc:

************************************************
Creating new SunPy Subclasses (Maps, TimeSeries)
************************************************

Writing a new Instrument Map Class
==================================

All instrument Map classes are subclasses of the generic `~sunpy.map.GenericMap` subclass.
`~sunpy.map.GenericMap` expect to be provided the data array as well as header dictionary and will parse the header metadata to the best of its ability based on common standards.
The instrument subclass implements the instrument-specific code to parse the metadata, apply any necessary procedures on the data array, as well as defining other things such what color tables to use.

In practice, the instrument subclass is not directly accessed by users.
The `~sunpy.map.Map` factory is the primary interface for creating Map objects.
Any subclass of `~sunpy.map.GenericMap` which defines a method named `~sunpy.map.GenericMap.is_datasource_for` will automatically be registered with the `~sunpy.map.Map` factory.
The ``is_datasource_for`` method is used by the `~sunpy.map.Map` factory to check if a file should use a particular instrument Map class.
This function can run any test it needs to determine this.
For example, it might check the value of the ``INSTRUMENT`` key in the metadata dictionary.
The following example shows how this works and includes a sample doc string that is compatible with the :ref:`Docs Guidelines for Data Sources`.

.. code-block:: python

    import sunpy.map
    class NextGenerationTelescopeMap(sunpy.map.GenericMap):
      "NextGenerationTelescope Map.

      The Next Generation Telescope is a type A telescope on board the XYZ mission.
      It operates in low Earth orbit with an altitude of 600 kmn and an inclination of 28.5 degrees.
      It is designed to observe the aliens on the Sun that are responsible for triggering the impulsive release of magnetic energy in the solar corona.
      It observes in the following 3 different passband in visible light, wavelength A, wavelength B, wavelength C.
      The primary emission processes in these passbands are process A and process B.

      The focal plane consists of a MAGIC detector with 2 x 2 pixels.
      The plate scale is 500 arcsec per pixel.
      The field of view is the whole Sun (1000 x 1000 arsec).
      It makes images in each passband every 10 minutes except for when it is in eclipse which occurs every approximately 30 minutes.

      It began operating on 2100 November 1.

      Notes
      -----
      Due to rise of our new insect overlords, the telescope was not operational from 2200 Jan to 2300 Jan.

      References
      ----------
      * List of all required references
      "

        def __init__(self, data, header, **kwargs):

            # will process the header according to common standards
            super(FutureMap, self).__init__(data, header, **kwargs)

            # Any NextGenerationTelescope Instrument-specific manipulation

        # used by the Map factory to determine if this subclass should be used
        @classmethod
        def is_datasource_for(cls, data, header, **kwargs):
            """Determines if data, header corresponds to a NextGenerationTelescope image"""
            # returns True only if this is data and header from NextGenerationTelescope
            return header.get('instrume', '').startswith('NextGenerationTelescope')


Writing a new Instrument TimeSeries Class
=========================================

To be written.
