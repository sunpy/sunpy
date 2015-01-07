.. _map:

=========
SunPy map
=========

.. module:: sunpy.map

.. testsetup::

    import sunpy

.. currentmodule:: sunpy.map

Overview
--------
One of core classes in SunPy is a Map. A SunPy Map object is simply a 
spatially-aware data array, often an image. In order to make it easy to work
with image data in SunPy, the Map object provides a number of methods for
commonly performed operations.

2D map objects are subclasses of sunpy.map.MapBase and all Map objects are 
created using the Map factory sunpy.map.Map.

.. Todo:
    1. Map factory and registration
    2. MapBase and Generic Map
    3. MapMeta and the seperation from the file io


Creating Map Objects
--------------------
SunPy Map objects are constructed using the special factory 
class :class:`Map`: ::

>>> x = sunpy.map.Map('file.fits')

The result of a call to `Map` will be either a `mapbase.GenericMap` object, 
or a subclass of `mapbase.GenericMap` which either deals with a specific type of data, 
e.g. `AIAMap` or `LASCOMap`, or if no instrument matches, a 2D map `mapbase.GenericMap`.


.. autoclass:: sunpy.map.Map
   
Map Classes
-----------
There are a series of base map classes which are specalised for each 
instrument. These subclass GenericMap and then register with
the Map factory class, which will direct instatiation of an instrument class if the correct
parameters are met. 

.. automodapi:: sunpy.map
    :no-main-docstr:
    :no-heading:

Writing a new Map Class
-----------------------

Map classes can be registered with the Map factory, even if the new class is not
officially part of SunPy.  This is good for prototyping new instruments.  For
example, to add a Map type for a future instrument, consider the code skeleton ::

    import sunpy
    
    class FutureMap(sunpy.GenericMap):
    
        def __init__(self, data, header, **kwargs):
        
            GenericMap.__init__(self, data, header, **kwargs)
            
            # Any Future Instrument specific keyword manipulation
       
       # Specify a classmethod that determines if the data-header pair matches
       # the new instrument
       @classmethod
       def is_datasource_for(cls, data, header, **kwargs):
            """Determines if header corresponds to an AIA image"""
            return header.get('instrume', '').startswith('FUTURESCOPE')
            
Then, to be able to instantiate a FutureMap using the Map() factory, one must
register the FutureMap type with the factory ::

    sunpy.map.Map.register(FutureMap, FutureMap.is_datasource_for)
    
If this line is placed correctly, for example in your subpackages __init__.py,
it can be guaranteed that the FutureMap is always accessible when your package
is imported.
