.. _sunpy-how-to-read-an-asdf-file:

****************************
Read an ASDF file into a Map
****************************

`ASDF <https://asdf-standard.readthedocs.io/en/latest/>`__ is a modern file format designed to meet the needs of the astronomy community [citation needed].
It has deep integration with Python, sunpy, and Astropy, as well as implementations in other languages.
It can be used to store known Python objects in a portable, well defined file format.
It is primarily useful for storing complex Astropy and sunpy objects in a way that can be loaded back into the same form as they were saved.
It is designed to be an archive file format, with human readable metadata and a simple on-disk layout.

sunpy currently implements support for saving `Map <sunpy.map.GenericMap>` and `coordinate frame <sunpy.coordinates.frames>` objects into asdf files.
As asdf tightly integrates into Python, saving a map to an asdf file will save the metadata, data and mask.
In comparison, the mask is not currently saved to FITS.

The following code shows to to save and load a sunpy Map to an asdf file:

.. code-block:: python

    >>> import asdf

    >>> import sunpy.map
    >>> from sunpy.data.sample import AIA_171_IMAGE  # doctest: +REMOTE_DATA

    >>> aiamap = sunpy.map.Map(AIA_171_IMAGE)  # doctest: +REMOTE_DATA
    >>> tree = {'amap': aiamap}  # doctest: +REMOTE_DATA
    >>> with asdf.AsdfFile(tree) as asdf_file:  # doctest: +REMOTE_DATA
    ...     asdf_file.write_to("sunpy_map.asdf")  # doctest: +REMOTE_DATA
    >>> input_asdf = asdf.open("sunpy_map.asdf")  # doctest: +REMOTE_DATA
    >>> input_asdf['amap']  # doctest: +REMOTE_DATA
        <sunpy.map.sources.sdo.AIAMap object at ...>
        SunPy Map
        ---------
        Observatory:                 SDO
        Instrument:          AIA 3
        Detector:            AIA
        Measurement:                 171.0 Angstrom
        Wavelength:          171.0 Angstrom
        Observation Date:    2011-06-07 06:33:02
        Exposure Time:               0.234256 s
        Dimension:           [1024. 1024.] pix
        Coordinate System:   helioprojective
        Scale:                       [2.402792 2.402792] arcsec / pix
        Reference Pixel:     [511.5 511.5] pix
        Reference Coord:     [3.22309951 1.38578135] arcsec
        array([[ -95.92475  ,    7.076416 ,   -1.9656711, ..., -127.96519  ,
                -127.96519  , -127.96519  ],
            [ -96.97533  ,   -5.1167884,    0.       , ...,  -98.924576 ,
                -104.04137  , -127.919716 ],
            [ -93.99607  ,    1.0189276,   -4.0757103, ...,   -5.094638 ,
                -37.95505  , -127.87541  ],
            ...,
            [-128.01454  , -128.01454  , -128.01454  , ..., -128.01454  ,
                -128.01454  , -128.01454  ],
            [-127.899666 , -127.899666 , -127.899666 , ..., -127.899666 ,
                -127.899666 , -127.899666 ],
            [-128.03072  , -128.03072  , -128.03072  , ..., -128.03072  ,
                -128.03072  , -128.03072  ]], dtype=float32)
    >>> input_asdf.close()  # doctest: +REMOTE_DATA

When saving a Map to ASDF, all maps are saved as a `.GenericMap` and not a specific source class.
This comes with some trade-offs.
If you are using custom map sources defined outside of the `sunpy` core package, and these sources are imported after asdf has been invoked for the first time (used, not just imported), then they will not be registered with the asdf converter.
Also if the custom map subclass is not registered with `sunpy.map.Map` upon loading of the map, it will be returned as a `.GenericMap`.
This approach has been chosen despite these limitations so that once a map is saved to an ASDF file it can always be loaded back into a map rather than the asdf library returning it as a Python dictionary.
It also follows the philosophy of the way maps are saved and loaded in the FITS format, where the components of the Map are serialized and the way meta data is handled depends solely on the contents of the ``.meta`` attribute.
