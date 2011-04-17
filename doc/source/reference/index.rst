===============
SunPy Reference
===============
Currently SunPy development is broken up between work on foundational classes
and routines such as the :class:`sunpy.Map` class, and a lose collection of
experimental code in the :mod:`sunpy.dev` module.

Classes
-------

Map
^^^
Maps are created in SunPy using the `Map()` function, which accepts either a
filename or a data array as input and returns either a generic Map object 
(BaseMap) or a subclass of BaseMap which deals with a specific type of data,
e.g. "AIAMap" or "LASCOMap".

.. automodule:: sunpy.data.map
.. automodule:: sunpy.data.BaseMap
   :members:
   :show-inheritance:

MapCube
^^^^^^^
MapCubes are similar to Map object except that they contain multiple 2d data
arrays, referred to as 'slices.'

.. automodule:: sunpy.data.MapCube
   :members:
   :show-inheritance:
   
Sun
^^^

.. automodule:: sunpy.Sun
   :members:
   :show-inheritance:
   
Experimental
------------
The purpose of the dev module is to provide a sandboxed environment where SunPy
developers can try out new code without having to decide on a formal interface.
Functions here may be incorporated into SunPy in the future as part of
the main namespace. Once SunPy is further along in it's development it is likely
that this module will be dropped altogether.

.. warning::

 Users should not write code for the main SunPy namespace that depends on
 functionality in the dev module. There is no guarantee that functions in the
 :mod:`sunpy.dev` module will be available in future versions of SunPy.

.. automodule:: sunpy.dev
    :members:
    
cm
^^
Experimental colormap functions

.. automodule:: sunpy.dev.cm
    :members:
