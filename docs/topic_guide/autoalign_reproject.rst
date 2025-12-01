.. _autoalign_vs_reproject:

*****************************************
Autoalignment vs. Reprojection
*****************************************

When working with solar data, you will often need to display or analyze maps in coordinate frames different from the one in which the data was observed. It is crucial to understand the difference between **autoalignment** (a visualization feature) and **reprojection** (a data manipulation technique).

Autoalignment
=============

Autoalignment is strictly a **visualization** technique. When you call :meth:`sunpy.map.GenericMap.plot` with a WCS that differs from the map's native coordinates, SunPy "autoaligns" the image by transforming the pixel coordinates to the display frame.

* **What it does:** It visually warps the image on your screen (e.g., by drawing each pixel as a quadrilateral) to match the axes.
* **What it does NOT do:** It does **not** modify the underlying ``.data`` array of the Map.
* **Performance:** Fast for interactive plots because no interpolation is required.
* **Default Behavior:** ``map.plot()`` automatically detects if autoalignment is necessary and applies it by default.

Use autoalignment when you simply want to **see** the data in a different orientation (e.g., looking at a rotated Sun) but do not need to perform calculations on the rotated data.

Reprojection
============

Reprojection is a **data transformation** technique. It uses the :meth:`sunpy.map.GenericMap.reproject_to` method to mathematically resample the map onto a new grid of pixels.

* **What it does:** It creates a brand new ``Map`` object. It calculates the value of every new pixel based on interpolation of the original data.
* **The Result:** The ``.data`` array of the new map is fundamentally different from the original.
* **Use Case:** Use reprojection when you need to perform pixel-by-pixel arithmetic (e.g., adding two maps together) or stack data that was taken from different viewpoints.

Visual Comparison
=================

The following example illustrates the difference. We take a map and rotate the view by 45 degrees.
Notice that **Autoalignment** keeps the data unchanged (visual rotation only), while **Reprojection** creates a new grid (notice the black padded corners where data was interpolated).

.. plot::
    :include-source:

    import matplotlib.pyplot as plt
    import astropy.units as u
    from astropy.wcs import WCS
    from sunpy.map import Map
    from sunpy.data.sample import AIA_171_IMAGE
    from sunpy.map.header_helper import make_fitswcs_header

    # Load the sample map
    m = Map(AIA_171_IMAGE)

    # Create a new Header that is rotated by 45 degrees
    target_header = make_fitswcs_header(
        m.data.shape,
        m.center,
        scale=u.Quantity(m.scale),
        rotation_angle=45 * u.deg,
        instrument="Rotated View"
    )
    target_wcs = WCS(target_header)

    fig = plt.figure(figsize=(12, 5))

    # Plot 1: Autoalignment
    ax1 = fig.add_subplot(1, 2, 1, projection=target_wcs)
    m.plot(axes=ax1, title="Autoalignment\n(Data Unchanged)")
    m.draw_grid(axes=ax1, color='w', alpha=0.5)

    # Plot 2: Reprojection
    m_reprojected = m.reproject_to(target_header)
    ax2 = fig.add_subplot(1, 2, 2, projection=m_reprojected)
    m_reprojected.plot(axes=ax2, title="Reprojection\n(New Data Array)")
    m_reprojected.draw_grid(axes=ax2, color='w', alpha=0.5)

    plt.show()

Summary
=======

.. list-table::
   :widths: 25 25 50
   :header-rows: 1

   * - Feature
     - Modifies ``.data``?
     - Best Used For
   * - **Autoalignment**
     - No
     - Visualization, Quick inspection, Overlays
   * - **Reprojection**
     - Yes (New Array)
     - Data Analysis, Map Arithmetic, Stacking
