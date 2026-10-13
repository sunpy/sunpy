.. _sunpy-topic-guide-create-new-map-class:

**********************************
Writing a new instrument Map class
**********************************

All instrument Map classes are subclasses of the generic `~sunpy.map.GenericMap` subclass.
`~sunpy.map.GenericMap` expects to be provided the data array and a header dictionary and will parse the header metadata to the best of its ability based on common standards.
The instrument subclass implements the instrument-specific code to parse the metadata, apply any necessary procedures on the data array, as well as defining other things such what color tables to use.

In practice, the instrument subclass is not directly accessed by users.
The `~sunpy.map.Map` factory is the primary interface for creating Map objects.
Any subclass of `~sunpy.map.GenericMap` which defines a method named ``is_datasource_for()`` will automatically be registered with the `~sunpy.map.Map` factory.
The ``is_datasource_for()`` method is used by the `~sunpy.map.Map` factory to check if a file should use a particular instrument Map class.
This function can run any test it needs to determine this.
For example, it might check the value of the "INSTRUMENT" key in the metadata dictionary.
The following example shows how this works and includes a sample doc string that is aligned with `the documentation guidelines <https://docs.sunpy.org/en/latest/dev_guide/contents/documentation.html#docs-guidelines-for-data-sources>`__:

.. code-block:: python

    import sunpy.map

    class NextGenerationTelescopeMap(sunpy.map.GenericMap):
      """
      NextGenerationTelescope Map.

      The Next Generation Telescope is a optical telescope on board the new space mission.
      It operates in low Earth orbit with an altitude of 600 km and an inclination of 28.5 degrees.
      It is designed to observe the mechanisms that are responsible for triggering the impulsive release of magnetic energy in the solar corona.
      It observes in the following 3 different passbands, in visible light, wavelength A, wavelength B, wavelength C.

      The focal plane consists of a detector with 4k x 4k pixels.
      The plate scale is 0.1 arcsec per pixel.
      The field of view is the whole Sun (1000 x 1000 arcsec).
      It makes images in each passband every 1 second except for when it is in eclipse which occurs every approximately 80 minutes.

      It began operating on 2100 November 1.

      Notes
      -----
      Due to failure of the filter wheel, 2 of the different passbands are no longer functional.

      References
      ----------
      * Mission Paper
      * Instrument Paper(s)
      * Data Archive
      * Mission documents
      """

        def __init__(self, data, header, **kwargs):

            # Will process the header according to FITS common standards
            super().__init__(data, header, **kwargs)

            # Any NextGenerationTelescope Instrument-specific manipulation.
            # Any metadata changes should be done by overloading
            # the corresponding attributes/methods.

        # Used by the Map factory to determine if this subclass should be used
        @classmethod
        def is_datasource_for(cls, data, header, **kwargs):
            """
            Determines if data, header corresponds to a NextGenerationTelescope image
            """
            # Returns True only if this is data and header from NextGenerationTelescope
            return header.get('instrume', '').startswith('NextGenerationTelescope')
