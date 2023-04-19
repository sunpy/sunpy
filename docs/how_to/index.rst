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
     - `my_map = Map('file.fits') <sunpy.map>`
   * - save a map to a FITS file
     - `my_map.save('another_file.fits') <sunpy.map.GenericMap.save>`
   * - get a quicklook summary of a map
     - `my_map.quicklook() <sunpy.map.GenericMap.quicklook>`
   * - plot a map
     - `my_map.plot() <sunpy.map.GenericMap.plot>`
   * - access the underlying data array of a map
     - `my_map.data <sunpy.map.GenericMap.data>`
   * - access the map metadata
     - `my_map.meta <sunpy.map.GenericMap.meta>`
   * - access the observer location
     - `my_map.observer_coordinate <sunpy.map.GenericMap.observer_coordinate>`
   * - remove the roll angle
     - `my_map.rotate() <sunpy.map.GenericMap.rotate>`
   * - plot the solar limb
     - `my_map.draw_limb() <sunpy.map.GenericMap.draw_limb>`
   * - overlay a heliographic grid
     - `my_map.draw_grid() <sunpy.map.GenericMap.draw_grid>`
   * - create a time series
     - `my_timeseries = TimeSeries('file.nc') <sunpy.timeseries>`
   * - plot a time series
     - `my_timeseries.plot() <sunpy.timeseries.GenericTimeSeries.plot>`
   * - concatenate two time series together
     - `my_timeseries.concatenate(another_timeseries) <sunpy.timeseries.GenericTimeSeries.concatenate>`
   * - parse the string representation of a timestamp
     - :func:`~sunpy.time.parse_time`
   * - calculate the Carrington rotation number at a given time
     - :func:`~sunpy.coordinates.sun.carrington_rotation_number`
   * - calculate the time corresponding to a given Carrington rotation
     - :func:`~sunpy.coordinates.sun.carrington_rotation_time`
