Experimental
============
The purpose of the dev module is to provide a sandboxed environment where SunPy
developers can try out new code without having to decide on a formal interface.
Functions here may be incorporated into SunPy in the future as part of
the main namespace. Once SunPy is further along in it's development it is likely
that this module will be dropped altogether.

.. warning::

 Users should not write code for the main SunPy namespace that depends on
 functionality in the  :mod:`sunpy.dev` module. There is no guarantee that 
 functions in the sunpy.dev module will be available in future versions of SunPy.

.. automodule:: sunpy.dev
    :members:
    
cm
--
Experimental colormap functions

.. automodule:: sunpy.dev.cm
    :members: