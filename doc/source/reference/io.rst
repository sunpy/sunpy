.. _io:

============
Input/Output
============

The I/O mdule provides file handling routines for sunpy.
The purpose of the I/O module is to handle all file methods so all other sunpy 
code can call sunpy.io. I/O implements file agnostic functions in io.
file_tools for reading a writing files, that automatically detect the filetype 
provided. Each file type supported by sunpy has a submodule in io which must 
implement read, get_header and write as a minimum. For each file type 
get_header should return an instance of sunpy.io.header.FileHeader which is a 
subclass of OrderedDict the standard base for all meta data in sunpy.

.. currentmodule:: sunpy

.. autosummary::
   :toctree: generated/

   io
   io.fits
   io.jp2