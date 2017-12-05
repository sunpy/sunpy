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

    >>> x = sunpy.map.Map('file.fits')  # doctest: +SKIP

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
Defined in `sunpy.map.sources` are a set of `~sunpy.map.GenericMap` subclasses
which convert the specific metadata and other differences in each instruments
data to the standard `~sunpy.map.GenericMap` interface.
These 'sources' also define things like the colormap and default
normalisation for each instrument.
These subclasses also provide a method, which describes to the `Map <sunpy.map.map_factory.MapFactory>` factory
which data and metadata pairs match its instrument.

.. automodapi:: sunpy.map
    :no-main-docstr:
    :no-heading:

.. _map-sources:

Instrument Map Classes
----------------------

.. automodapi:: sunpy.map.sources


Writing a new Instrument Map Class
----------------------------------

Any subclass of `~sunpy.map.GenericMap` which defines a method named
`~sunpy.map.GenericMap.is_datasource_for` will automatically be registered with
the `Map <sunpy.map.map_factory.MapFactory>` factory. The ``is_datasource_for`` method describes the form of the
data and metadata for which the `~sunpy.map.GenericMap` subclass is valid. For
example it might check the value of the ``INSTRUMENT`` key in the metadata
dictionary.
This makes it straightforward to define your own
`~sunpy.map.GenericMap` subclass for a new instrument or a custom data source
like simulated data. These classes only have to be imported for this to work, as
demonstrated by the following example.

.. code-block:: python

    import sunpy.map
    class FutureMap(sunpy.map.GenericMap):

        def __init__(self, data, header, **kwargs):

            super(FutureMap, self).__init__(data, header, **kwargs)

            # Any Future Instrument specific keyword manipulation

       # Specify a classmethod that determines if the data-header pair matches
       # the new instrument
       @classmethod
       def is_datasource_for(cls, data, header, **kwargs):
            """Determines if header corresponds to an AIA image"""
            return header.get('instrume', '').startswith('FUTURESCOPE')


This class will now be available through the `Map <sunpy.map.map_factory.MapFactory>` factory as long as this
class has been defined, i.e. imported into the current session.

If you do not want to create a method named ``is_datasource_for`` you can
manually register your class and matching method using the following method

.. code-block:: python

    import sunpy.map

    sunpy.map.Map.register(FutureMap, FutureMap.some_matching_method)
