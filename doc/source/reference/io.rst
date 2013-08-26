.. _io:

============
Input/Output
============

The I/O module provides file handling routines for sunpy.
The purpose of the I/O module is to handle all file methods so all other sunpy 
code can call sunpy.io. I/O implements file agnostic functions in io.file_tools 
for reading a writing files, that automatically detect the filetype 
provided. Each file type supported by sunpy has a submodule in io which must 
implement `read`, `get_header` and `write` as a minimum. 

A file in SunPy is defined in two parts, data and meta data. Each implemented 
file type should provide support for both. get_header should return the meta 
data which should be a instance of sunpy.io.header.FileHeader or a subclass 
there of if any conversion is required to validate a header. *i.e.* truncating 
keys in FITS files. sunpy.io.header.FileHeader is a subclass of OrderedDict 
which is the standard base for all SunPy meta data classes.

To enable support for file types that include multiple data, header 
combinations `read` and `get_header` must both return lists even if there is 
only one item in the list.

The main functions avalible through the sunpy.io namespace are as follows:
 
.. automodule:: sunpy.io.file_tools
   :members:

The I/O submodule currently supports two file types FITS and JPEG2000.
The documentation for these submodules is avalible in the table below.

.. currentmodule:: sunpy

.. autosummary::
   :toctree: generated/

   io.fits
   io.jp2
   io.header
