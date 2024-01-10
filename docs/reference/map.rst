.. _map:

Maps (`sunpy.map`)
******************

.. module:: sunpy.map

sunpy Map objects are constructed using the special factory class:

.. autoclass:: sunpy.map.map_factory.MapFactory

.. _map-classes:

.. currentmodule:: sunpy.map

sunpy.map Package
=================

Map objects are constructed using the special factory class: `~sunpy.map.Map`.
All sunpy Maps are derived from `sunpy.map.GenericMap`, all the methods and attributes are documented in that class.

The result of a call to `~sunpy.map.Map` will be either a `~sunpy.map.mapbase.GenericMap` object if no instrument matches, or a subclass of `~sunpy.map.mapbase.GenericMap` which deals with a specific source of data, e.g., `~sunpy.map.sources.sdo.AIAMap` or `~sunpy.map.sources.soho.LASCOMap` (see :ref:`map-classes` to see a list of all of them).

.. automodapi:: sunpy.map
    :no-main-docstr:
    :inherited-members:
    :include-all-objects:
    :skip: make_fitswcs_header
    :skip: get_observer_meta
    :skip: make_heliographic_header

Header helpers
==============

The header_helper sub-module contains helper functions for generating FITS-WCS headers from Python objects.

.. automodapi:: sunpy.map.header_helper
    :no-main-docstr:
    :inherited-members:
    :include-all-objects:

.. _map-sources:

Instrument Map Classes
======================
Defined in ``sunpy.map.sources`` are a set of `~sunpy.map.GenericMap` subclasses which convert the specific metadata and other differences in each instruments data to the standard `~sunpy.map.GenericMap` interface.
These 'sources' also define things like the colormap and default normalization for each instrument.
These subclasses also provide a method, which describes to the `~sunpy.map.Map` factory which data and metadata pairs match its instrument.

.. automodapi:: sunpy.map.sources
    :no-main-docstr:
    :inherited-members:
