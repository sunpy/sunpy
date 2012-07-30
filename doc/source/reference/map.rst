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
commonly performed operations. Further, because SunPy Map objects are
instances inherit from the NumPy `ndarray <http://docs.scipy.org/doc/numpy/reference/generated/numpy.ndarray.html>`_ 
datatype, they behave like ndarrays and support the same operations as ndarrays.

Creating Map Objects
^^^^^^^^^^^^^^^^^^^^
SunPy Map objects are constructed using the special function :func:`make_map`: ::

>>> x = sunpy.make_map('file.fits')

The result of a call to `make_map` will be either a generic `Map` object, 
or a subclass of `Map` which deals with a specific type of data, e.g. 
`AIAMap` or `LASCOMap`.

.. autofunction:: make_map
   
Map Classes
^^^^^^^^^^^

:class:`Map`
""""""""""""""""
The top-level class from which all other Maps inherit from.

.. autoclass:: Map

    
:class:`MapHeader`
""""""""""""""""""

Header information for `Map` objects are stored in a class called 
:class:`MapHeader`.

.. autoclass:: MapHeader

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

