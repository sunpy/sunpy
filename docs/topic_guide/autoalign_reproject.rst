.. _sunpy-topic-guide-map-autoalign-reproject:

*******************************************
Autoalignment vs. Reprojection in Map Plots
*******************************************

When displaying a `~sunpy.map.GenericMap` on coordinate axes that differ from the map's native coordinate system, each map may have its own coordinate frame, projection, and pixel grid.
To align them visually or for analysis, there are two fundamentally different approaches: **autoalignment** and **reprojection**.
This page explains the difference, when to use each, and how they affect your data.

Quick summary
=============

Autoalignment is a visualization technique that does not modify the map's pixel data.
Reprojection is a data transformation that produces a new map by resampling data onto a new pixel grid.

Autoalignment
=============

Autoalignment is controlled by the ``autoalign`` parameter of :meth:`~sunpy.map.GenericMap.plot`.
When the WCS of the map differs from the WCS of the `~astropy.visualization.wcsaxes.WCSAxes` axes, autoalignment draws the map's pixels in a coordinate-aware fashion so that they appear correctly aligned on the axes.
Crucially, the ``.data`` array of the map is **never modified**.

There are three modes:

* ``autoalign='mesh'`` draws each map pixel individually as a quadrilateral.
  This is the most general approach and works even when the coordinate transformation produces non-convex or highly warped results.
* ``autoalign='image'`` draws the entire map as a single warped image using an affine transform.
  This is usually faster than the mesh-based approach, but has limitations when the transformation is highly non-linear.
* ``autoalign=True`` (the default since sunpy 7.0) automatically determines whether to use the mesh-based or image-based approach depending on the coordinate transformation involved.

.. note::

    Autoalignment is purely a visualization tool.
    It does not change the underlying ``.data`` array or the WCS of the map.
    If you need the pixel data itself to be on a new grid — for example, to compute a difference image or co-add observations — use reprojection instead.

.. note::

    When plotting a `~sunpy.coordinates.Helioprojective` map with autoalignment, off-disk coordinates may not be visible because they are undefined without a screen assumption.
    See `~sunpy.coordinates.SphericalScreen` for how to handle off-disk data in this context.

Reprojection
============

Reprojection transforms the map's pixel data onto a new pixel grid defined by a target WCS.
This is done using :meth:`~sunpy.map.GenericMap.reproject_to`, which returns a **new** `~sunpy.map.GenericMap` with a resampled ``.data`` array.
The original map is not modified.

Unlike autoalignment, reprojection involves interpolation: each pixel in the output grid is computed by sampling the input data at the corresponding world coordinate.
This means the output data array is genuinely different from the input, and the two maps share an identical pixel grid after reprojection.

Since sunpy 7.0, :meth:`~sunpy.map.GenericMap.reproject_to` supports automatic determination of the output extent via the ``auto_extent`` keyword argument.
Since sunpy 7.1, it also accepts a ``preserve_date_obs`` keyword argument to retain the original observation time of the map in the reprojected result.
This is useful when reprojecting a series of images to a common coordinate frame while preserving each image's timestamp.

Use reprojection when you need to:

* Perform pixel-level arithmetic between two maps (e.g., computing a ratio or difference image).
* Co-add or stack images from different viewpoints or projections.
* Ensure that two maps share an identical pixel grid before analysis.

When to use which
=================

.. list-table::
   :widths: 20 20 60
   :header-rows: 1

   * - Approach
     - Modifies ``.data``?
     - Use when...
   * - Autoalignment
     - No
     - You want to **display** the map in a different coordinate frame or projection without changing the data.
   * - Reprojection
     - Yes (new map)
     - You need the pixel data itself to be on a new grid, e.g., for arithmetic, co-adding, or stacking.

See also
========

* :ref:`sunpy-tutorial-maps` — introductory tutorial covering `~sunpy.map.GenericMap` basics
* :ref:`sunpy-tutorial-map-plotting-maps` — tutorial section on visualizing maps
* :ref:`sunpy-how-to-index` — short task-focused guides

Examples
========

.. minigallery:: ../examples/map_transformations/autoalign_aia_hmi.py ../examples/map_transformations/reprojection_align_aia_hmi.py
