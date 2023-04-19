.. _how_to_guide:

*************
How-To Guides
*************

These guides provide examples of how to perform specific tasks with sunpy.
These are recipes that do not provide much in depth explanation and assume you have some knowledge of what ``sunpy`` is and how it works.
If you're starting fresh you might want to check out the :ref:`tutorial` first.

.. toctree::
   :maxdepth: 1

   search_multiple_wavelengths.rst


Quick Reference
---------------

The following table is meant as a quick reference.
For more complete code examples, see the how-to guides above.

.. list-table::
   :widths: auto
   :header-rows: 1

   * - How do I...
     - Solution
   * - create a map from a FITS file
     - :mod:`Map('file.fits') <sunpy.map>`
   * - save a map to a FITS file
     - `Map.save() <sunpy.map.GenericMap.save>`
   * - get a quicklook summary of a map
     - `Map.quicklook() <sunpy.map.GenericMap.quicklook>`
   * - plot a map
     - `Map.plot() <sunpy.map.GenericMap.plot>`
   * - access the underlying data array of a map
     - `Map.data <sunpy.map.GenericMap.data>`
   * - access the metadata of a map
     - `Map.meta <sunpy.map.GenericMap.meta>`
   * - access the observer location of a map
     - `Map.observer_coordinate <sunpy.map.GenericMap.observer_coordinate>`
   * - find the pixel in a map corresponding to a world coordinate
     - `Map.world_to_pixel() <sunpy.map.GenericMap.world_to_pixel>`
   * - derotate a map to remove the roll angle
     - `Map.rotate() <sunpy.map.GenericMap.rotate>`
   * - resample a map to a different resolution
     - `Map.resample() <sunpy.map.GenericMap.resample>`
   * - reproject a map to a different coordinate system
     - `Map.reproject_to() <sunpy.map.GenericMap.reproject_to>`
   * - plot the solar limb on top of a map
     - `Map.draw_limb() <sunpy.map.GenericMap.draw_limb>`
   * - overlay a heliographic grid on top a map
     - `Map.draw_grid() <sunpy.map.GenericMap.draw_grid>`
   * - plot a time series
     - `TimeSeries.plot() <sunpy.timeseries.GenericTimeSeries.plot>`
   * - concatenate two time series together
     - `TimeSeries.concatenate() <sunpy.timeseries.GenericTimeSeries.concatenate>`
   * - parse the string representation of a timestamp
     - :func:`~sunpy.time.parse_time`
   * - calculate the Carrington rotation number at a given time
     - :func:`~sunpy.coordinates.sun.carrington_rotation_number`
   * - calculate the time corresponding to a given Carrington rotation
     - :func:`~sunpy.coordinates.sun.carrington_rotation_time`
