.. _map:

===
Map
===

.. currentmodule:: sunpy

Creating Map Objects
^^^^^^^^^^^^^^^^^^^^

SunPy Map objects are constructed using the special function :func:`Map`: ::

>>> x = sunpy.Map('file.fits')

The result of a call to `Map` will be either a generic `BaseMap` object, or a 
subclass of `BaseMap` which deals with a specific type of data, e.g. 
`AIAMap` or `LASCOMap`.

.. autosummary::
   :toctree: generated/

   Map
   
BaseMap Class
^^^^^^^^^^^^^

The top-level class from which all other Maps inherit from is `BaseMap`.

.. autosummary::
   :toctree: generated/

    map.BaseMap
    
MapHeader
^^^^^^^^^

Header information for `Map` objects are stored in a class called 
:class:`MapHeader`.

.. autosummary::
   :toctree: generated/

    MapHeader
