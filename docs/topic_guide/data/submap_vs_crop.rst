.. _sunpy-topic-guide-submap-vs-crop:

*************************************
Cutting out a region: submap vs crop
*************************************

`~sunpy.map.GenericMap` provides two ways to cut out a region from a map; :meth:`~sunpy.map.GenericMap.submap` and the inherited :meth:`~ndcube.NDCube.crop`.
Both take world coordinates and return a rectangular block of pixels.
Given the same two coordinates, they often return the same region, but there are maps for which their returns will differ, for example when the pixel grid has been rotated relative to the world coordinate axes.
This guide explains how the methods interpret their inputs and therefore why those differences "crop" up.

What they have in common
========================

A map stores its data in a rectangular array; a cutout must contain a whole number of rows and columns.
Neither `submap` nor `crop` resamples or rotates the data, both convert the input coordinates to pixel positions and return the smallest rectangle of the map's pixels that contains them.

Where they differ
=================

The methods differ in how they use the supplied coordinates to define the pixel rectangle.

:meth:`~ndcube.NDCube.crop` uses only the coordinates you pass and makes no assumptions about their meaning.
You can pass any number of coordinates.

:meth:`~sunpy.map.GenericMap.submap` takes only two arguments and treats these as opposite corners of a rectangle *in world coordinates*, then derives the other two corners prior to cropping.
For example, when provided with world coordinates at ``(-900", -900")`` and ``(0", -200")`` this implies two additional corners: ``(-900", -200")`` and ``(0", -900")``.
It passes all four corners to :meth:`~ndcube.NDCube.crop`.

To get the same result with :meth:`~ndcube.NDCube.crop`, you could calculate those additional two corners and pass all four world-coordinate corners yourself.
:meth:`~sunpy.map.GenericMap.submap` simply derives the two additional corners for you.

For many maps, the pixel grid is aligned with the world coordinate axes: moving along a row changes longitude, and moving along a column changes latitude.
The rectangle defined by the two supplied corners is then aligned with the pixel grid. The derived corners fall on its edges, and both methods return identical regions.

On a rotated map
================

On a rotated map, the sides of a world-coordinate rectangle do not follow the pixel rows and columns, so a pixel-aligned cutout must extend beyond the world-coordinate rectangle's edges.

The two coordinates passed to :meth:`~ndcube.NDCube.crop` define one diagonal of the rectangle.
The rectangle produced by :meth:`~ndcube.NDCube.crop` is minimised to contain just those two points, while :meth:`~sunpy.map.GenericMap.submap` also includes the two derived corners, which extend the rectangle.
The resulting pixel bounding boxes are different:

.. code-block:: python

    >>> import astropy.units as u
    >>> from astropy.coordinates import SkyCoord
    >>> import sunpy.map
    >>> from sunpy.data.sample import SWAP_LEVEL1_IMAGE  # doctest: +SKIP

    >>> rotated = sunpy.map.Map(SWAP_LEVEL1_IMAGE).rotate(30 * u.deg)  # doctest: +SKIP
    >>> bottom_left = SkyCoord(-900 * u.arcsec, -900 * u.arcsec, frame=rotated.coordinate_frame)  # doctest: +SKIP
    >>> top_right = SkyCoord(0 * u.arcsec, -200 * u.arcsec, frame=rotated.coordinate_frame)  # doctest: +SKIP

    >>> rotated.crop(bottom_left, top_right).data.shape  # doctest: +SKIP
    (335, 137)
    >>> rotated.submap(bottom_left, top_right=top_right).data.shape  # doctest: +SKIP
    (335, 358)

The plot shows the tilted rectangle and the pixel bounding boxes:

.. plot::
    :context: close-figs

    import matplotlib.pyplot as plt
    import numpy as np
    from matplotlib.patches import Polygon, Rectangle

    import astropy.units as u
    from astropy.coordinates import SkyCoord

    import sunpy.data.sample
    import sunpy.map

    rotated = sunpy.map.Map(sunpy.data.sample.SWAP_LEVEL1_IMAGE).rotate(30 * u.deg)

    def to_pixel(lon, lat):
        coord = SkyCoord(lon * u.arcsec, lat * u.arcsec, frame=rotated.coordinate_frame)
        return np.array(rotated.wcs.world_to_pixel(coord))

    bl = to_pixel(-900, -900)
    tr = to_pixel(0, -200)
    br = to_pixel(0, -900)
    tl = to_pixel(-900, -200)

    def bounding_box(points, **kwargs):
        xs, ys = np.array(points).T
        return Rectangle((xs.min(), ys.min()), xs.max() - xs.min(), ys.max() - ys.min(),
                         facecolor="none", lw=2.5, **kwargs)

    fig, ax = plt.subplots(figsize=(6, 6))
    ax.add_patch(Polygon([bl, br, tr, tl], facecolor="0.93", edgecolor="0.5", lw=1, ls=":"))
    ax.add_patch(bounding_box([bl, tr, br, tl], edgecolor="tab:cyan", label="submap"))
    ax.add_patch(bounding_box([bl, tr], edgecolor="red", ls="--", label="crop(bl, tr)"))

    ax.plot(*np.array([bl, tr]).T, "o", color="red", ms=9, zorder=5)
    ax.plot(*np.array([br, tl]).T, "o", color="tab:cyan", mfc="white", ms=9, zorder=5)
    for name, point in [("bl", bl), ("tr", tr), ("br", br), ("tl", tl)]:
        ax.annotate(name, point, xytext=(8, 6), textcoords="offset points")

    ax.set_xlabel("pixel column")
    ax.set_ylabel("pixel row")
    ax.set_aspect("equal")
    ax.legend(loc="upper left")
    ax.margins(0.15)
    ax.grid(alpha=0.25)

    plt.show()

The grey polygon is the requested region. The red markers show the supplied ``bl`` and ``tr`` coordinates; the open markers show ``br`` and ``tl``, the two corners derived by :meth:`~sunpy.map.GenericMap.submap`.

The out from :meth:`~ndcube.NDCube.crop` is a narrow strip because the supplied diagonal runs steeply across the pixel grid.
The submap, whose bounds include the points from the other diagonal, extends the region horizontally.
Both cutouts have the same height because the supplied diagonal already spans the full height.

Which one to use
================

Use :meth:`~sunpy.map.GenericMap.submap` when you want to cut out the full rectangle on the sky defined by two opposite corners.
It includes that world-coordinate region regardless of the map's rotation, though the pixel cutout may extend beyond it.

``crop`` is more flexible, use :meth:`~ndcube.NDCube.crop` when you need to include a set of any length of coordinates, such as a series of features, without implying the other corners of a rectangle when using 2.
You can also always pass all four world-coordinate corners to define the rectangle explicitly with ``crop``.


Summary
=======

Given the same two coordinates, :meth:`~ndcube.NDCube.crop` uses a subset of the points used by :meth:`~sunpy.map.GenericMap.submap`.
The bounding box of a subset cannot be larger than that of the full set, so the region returned by :meth:`~ndcube.NDCube.crop` is always contained within the region returned by :meth:`~sunpy.map.GenericMap.submap`.

The methods return the same region when the pixel grid is aligned with the world coordinate axes, including after rotations by multiples of 90 degrees.
Otherwise, :meth:`~ndcube.NDCube.crop` returns a smaller region. This can happen with other rotations or with a skewed grid whose axes are not perpendicular.

This comparison applies only when both methods receive the same two coordinates.
With more points, :meth:`~ndcube.NDCube.crop` returns the smallest region containing them all, which may be larger.

See :ref:`sphx_glr_generated_gallery_map_creating_a_submap.py` and :ref:`sphx_glr_generated_gallery_map_cropping_a_map.py` for examples.
