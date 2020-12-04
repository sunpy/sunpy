SunPy IO
********

This submodule contains two types of routines, the first reads (data, header)
pairs from files in a way similar to FITS files. The other is special readers
for files that are commonly used in solar physics.

.. warning::

   When reading FITS files, it is strongly recommended that `astropy.io.fits` be used over the tools in `sunpy.io`.
   The sunpy FITS reader is designed to meet the needs of map, and does not represent the structure of the FITS file well.

Special File Readers
====================

.. _iospecialgenx:
.. automodapi:: sunpy.io.special.genx


asdf (Advanced Scientific Data Format)
--------------------------------------

`asdf <https://asdf.readthedocs.io/en/latest/>`__ is a modern file format
designed to meet the needs of the astronomy community. It has deep integration
with Python and SunPy and Astropy as well as implementations in other languages.
It can be used to store known Python objects in a portable, well defined file
format. It is primarily useful for storing complex Astropy and SunPy objects in
a way that can be loaded back into the same form as they were saved.

SunPy currently implements support for saving `Map <sunpy.map.GenericMap>` and
`coordinate frame <sunpy.coordinates.frames>` objects into asdf files. As asdf
tightly integrates into Python, saving a map to an asdf file will save the
metadata, data, mask and the shift. The mask and shift are not currently saved
to FITS. The following code shows to to save and load a SunPy Map to an asdf
file

.. doctest-requires:: asdf

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

Unified File Readers
====================

.. automodapi:: sunpy.io

.. automodapi:: sunpy.io.header

.. _iofits:
.. automodapi:: sunpy.io.fits

.. _iojp2:
.. automodapi:: sunpy.io.jp2

.. _ioana:
.. automodapi:: sunpy.io.ana
