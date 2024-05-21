Input/output (``sunpy.io``)
***************************

The primary focus of ``sunpy.io`` is to provide a common interface for reading and writing.
``sunpy.io`` contains readers for files that are commonly used in solar physics.

These include:

- Public API
  - GENX
  - ANA
  - NOAA SWPC Solar Region Summary (SRS)
  - ASDF

The other readers are intended for use by `sunpy.map` and `sunpy.timeseries` and are not intended to be used directly by users.

Special File Readers
====================

.. automodapi:: sunpy.io.special.genx

.. automodapi:: sunpy.io.ana

.. automodapi:: sunpy.io.special.srs

.. automodapi:: sunpy.io.special.asdf
