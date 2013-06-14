.. _map:

----------
SunPy Maps
----------

.. currentmodule:: sunpy.map

Overview
^^^^^^^^
One of core classes in SunPy is a Map. A SunPy Map object is simply a 
spatially-aware data array, often an image. In order to make it easy to work
with image data in SunPy, the Map object provides a number of methods for
commonly performed operations.

2D map objects are subclasses of sunpy.map.MapBase and all Map objects are 
created using the Map factory sunpy.Map.

Updated Map Layout
^^^^^^^^^^^^^^^^^^
Todo:

1. Map factory and registration
2. MapBase and Generic Map
3. MapMeta and the seperation from the file io


Creating Map Objects
^^^^^^^^^^^^^^^^^^^^
SunPy Map objects are constructed using the special factory 
class :class:`Map`: ::

>>> x = sunpy.Map('file.fits')

The result of a call to `Map` will be either a `MapBase` object, 
or a subclass of `MapBase` which either deals with a specific type of data, 
e.g. `AIAMap` or `LASCOMap`, or a 2D map `GenericMap`.

.. autoclass:: Map
   
Map Classes
^^^^^^^^^^^
There are a series of base map classes which are specalised for each 
instrument. These subclass one of the MapBase derivaties and then register with
 the Map factory class which will instancestate a instrument class if the 
parameters are met. 

:class:`MapBase`
""""""""""""""""
The top-level class from which all other ND Maps inherit from.

.. autoclass:: MapBase

:class:`GenericMap`
"""""""""""""""""""
This is the top level 2D map class, containg processing and visualisation 
routines designed to work with 2D data.

.. autoclass:: `GenericMap`
    
:class:`MapMeta`
""""""""""""""""""

Meta data for `MapBase` objects are stored in a class called 
:class:`MapMeta`.

.. autoclass:: MapMeta

:class:`CompositeMap`
"""""""""""""""""""""
A Composite Map is a Map object which contains one or more layers, representing
for example a stack of images with varying opacities.

.. autoclass:: CompositeMap

:class:`MapCube`
"""""""""""""""""""""
A MapCube is a three-dimension generalization of the Map class, for example,
a time series of images.

.. autoclass:: MapCube

