.. _sunpy-how-to-index:

*************
How-To Guides
*************

These how-to guides provide examples of how to perform specific tasks with sunpy.
They are recipes that do not provide much in-depth explanation and assume you have some knowledge of what sunpy is and how it works.
If you're starting fresh you might want to check out the :ref:`sunpy-tutorial-index` first.

.. toctree::
   :maxdepth: 1

   coord_components
   create_a_map
   create_coords
   create_custom_map
   create_custom_timeseries
   create_rectangle_on_map
   fix_map_metadata
   manipulate_grid_lines
   parse_time
   read_asdf_file
   remote_data_manager
   search_multiple_wavelengths
   search_vso
   transform_coords

**Quick Reference**

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
     - `my_map.quicklook() <sunpy.map.GenericMap.quicklook>`, `my_map.peek() <sunpy.map.GenericMap.peek>`
   * - plot a map
     - `my_map.plot() <sunpy.map.GenericMap.plot>`
   * - access the underlying data array of a map
     - `my_map.data <sunpy.map.GenericMap.data>`
   * - access the map metadata
     - `my_map.meta <sunpy.map.GenericMap.meta>`
   * - make a copy of the map data array
     - `my_map.data.copy() <numpy.ndarray.copy>`
   * - make a copy of the whole map
     - `copy.deepcopy(my_map) <copy.deepcopy>`
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
   * - convert a time series to a `pandas.DataFrame`
     - `my_timeseries.to_dataframe() <sunpy.timeseries.GenericTimeSeries.to_dataframe>`
   * - convert a time series to an `astropy.table.Table`
     - `my_timeseries.to_table() <sunpy.timeseries.GenericTimeSeries.to_table>`
   * - parse the string representation of a timestamp
     - :func:`~sunpy.time.parse_time`
   * - calculate the Carrington rotation number at a given time
     - :func:`~sunpy.coordinates.sun.carrington_rotation_number`
   * - calculate the time corresponding to a given Carrington rotation
     - :func:`~sunpy.coordinates.sun.carrington_rotation_time`
   * - see all of the available solar constants
     - :func:`sunpy.sun.constants.print_all`
