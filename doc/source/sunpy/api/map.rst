.. _map:

----
Maps
----

.. module:: sunpy.map

.. testsetup::

    import sunpy

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
e.g. `AIAMap` or `LASCOMap`, or if no instrument matches, a 2D map `GenericMap`.

The SunPy Map factory accepts a wide variety if inputs for creating maps::

* Preloaded tuples of (data, header) pairs

>>> mymap = sunpy.Map((data, header))

* data, header pairs, not in tuples

>>> mymap = sunpy.Map(data, header)

* File names

>>> mymap = sunpy.Map('file1.fits')

* All fits files in a directory by giving a directory

>>> mymap = sunpy.Map('local_dir/sub_dir')

* Some regex globs

>>> mymap = sunpy.Map('eit_*.fits')

* URLs

>>> mymap = sunpy.Map(url_str)

* Lists of any of the above

>>> mymap = sunpy.Map(['file1.fits', 'file2.fits', 'file3.fits', 'directory1/'])

* Any mixture of the above not in a list

>>> mymap = sunpy.Map((data, header), data2, header2, 'file1.fits', url_str, 'eit_*.fits')

.. method:: sunpy.map.Map.__new__

.. autoclass:: sunpy.map.Map
   
Map Classes
^^^^^^^^^^^
There are a series of base map classes which are specalised for each 
instrument. These subclass GenericMap and then register with
 the Map factory class, which will direct instatiation of an instrument class if the correct
parameters are met. 

:class:`GenericMap`
"""""""""""""""""""
This is the top level 2D map class, containing processing and visualisation 
routines designed to work with 2D data.

.. autoclass:: sunpy.map.GenericMap
    
:class:`MapMeta`
""""""""""""""""""

Meta data for `Map` objects are stored in a class called  :class:`MapMeta`.

.. autoclass:: sunpy.map.MapMeta

:class:`CompositeMap`
"""""""""""""""""""""
A Composite Map is a Map object which contains one or more layers, representing
for example a stack of images with varying opacities.

.. autoclass:: sunpy.map.CompositeMap

:class:`MapCube`
"""""""""""""""""""""
A MapCube is a three-dimension generalization of the Map class, for example,
a time series of images.

.. autoclass:: sunpy.map.MapCube



Writing a new Map Class
^^^^^^^^^^^^^^^^^^^^^^^

Map classes can be registered with the Map factory, even if the new class is not
officially part of SunPy.  This is good for prototyping new instruments.  For
example, to add a Map type for a future instrument, consider the code skeleton::

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

    sunpy.Map.register(FutureMap, FutureMap.is_datasource_for)
    
If this line is place correctly, for example in your subpackages __init__.py,
it can be guaranteed that the FutureMap is always accessible when your package
is imported.