.. _reference:

====================
SunPy Code Reference
====================

This document is the API reference for SunPy, it will eventually be complete. 
Currently the following submodules and functions are avalible from the top namespace:

Modules:
    * sunpy.map 
    * sunpy.sun 
    * sunpy.lightcurve 
Functions:
    * sunpy.map.Map 
    * sunpy.map.make_map 
    * sunpy.util.config.load_config 
    * sunpy.util.config.print_config 
All functions in: 
    * sunpy.cm 
sunpy.data: 
    * AIA_171_IMAGE 
    * RHESSI_IMAGE 
    * EIT_195_IMAGE 
    * RHESSI_EVENT_LIST 
    * CALLISTO_IMAGE 

Each submodule should be listed below and document every public facing
method.  Any change to this API should result in a major version
change.

.. module:: sunpy

.. toctree::
   :maxdepth: 2

   cm
   gui
   hek
   instr
   io
   lightcurve
   map
   rhessi
   sun
   time
   vso
   wcs
   spectrogram
   callisto
   util
   image
