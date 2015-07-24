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

2D map objects are subclasses of `~sunpy.map.MapBase` and all Map objects are
created using the Map factory `~sunpy.map.Map`.

A number of instrument are supported by subclassing this base object. See
:ref:`map-sources` to see a list of all of them. More complex subclasses are also
available. See :ref:`map-classes`.

.. todo :
    1. Map factory and registration
    2. MapBase and Generic Map
    3. MapMeta and the separation from the file io


Creating Map Objects
--------------------
SunPy Map objects are constructed using the special factory
class `~sunpy.map.Map`: ::

    x = sunpy.map.Map('file.fits')

The result of a call to `~sunpy.map.Map` will be either a `~sunpy.map.mapbase.GenericMap` object,
or a subclass of `~sunpy.map.mapbase.GenericMap` which either deals with a specific type of data,
e.g. `~sunpy.map.sources.sdo.AIAMap` or `~sunpy.map.sources.soho.LASCOMap`
(see :ref:`map-classes` to see a list of all of them), or if no
instrument matches, a 2D map `~sunpy.map.mapbase.GenericMap`.


.. autoclass:: sunpy.map.map_factory.MapFactory


Using Map Objects
-----------------

Once a map object has been created using `~sunpy.map.Map` it will be a instance
or a subclass of the `~sunpy.map.mapbase.GenericMap` class. Irrespective of
the instrument the map is constructed for, all maps behave the same and are
interchangeable with one another. It is possible to manipulate the map or access
meta data about the map from the methods and properties of the map class.
The following documentation of `~sunpy.map.mapbase.GenericMap` lists the
attributes and methods that are available on all Map objects.

.. autoclass:: sunpy.map.mapbase.GenericMap
   :members:

.. _map-classes:

Map Classes
-----------
There are a series of base map classes which are specialised for each
instrument. These subclass GenericMap and then register with
the Map factory class, which will be direct instantiation of an instrument class if the correct
parameters are met.

.. automodapi:: sunpy.map
    :no-main-docstr:
    :no-heading:

.. _map-sources:

Instrument Map Classes
----------------------

.. automodapi:: sunpy.map.sources


Writing a new Instrument Map Class
----------------------------------

Map classes can be registered with the Map factory, even if the new class is not
officially part of SunPy.  This is good for prototyping new instruments.  For
example, to add a Map type for a future instrument, consider this code skeleton::

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
register the FutureMap type with the factory::

    sunpy.map.Map.register(FutureMap, FutureMap.is_datasource_for)

If this line is placed correctly, for example in your subpackages `__init__.py`,
it can be guaranteed that the FutureMap is always accessible when your package
is imported.
