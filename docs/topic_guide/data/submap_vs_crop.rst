.. _sunpy-topic-guide-submap-vs-crop:

*************************************
Cutting out a region: submap vs crop
*************************************

`~sunpy.map.GenericMap` offers two ways to cut a smaller map from a larger one.
:meth:`~sunpy.map.GenericMap.submap` has long been part of SunPy, while :meth:`~ndcube.NDCube.crop` is inherited from `ndcube.NDCube`.
Both take world coordinates and return a rectangular block of pixels. On many maps, they return the same result when given the same two coordinates.
This guide explains how the methods differ and why their results can diverge.

What they have in common
========================

A map holds its data in a rectangular array; any cutout has to be a whole number of rows and columns.
Neither method resamples or rotates the data. Each converts the input coordinates to pixel positions and returns the smallest rectangle of the parent map's pixels that contains them.

Where they differ
=================

They differ in which coordinates they use to define the pixel rectangle.

:meth:`~ndcube.NDCube.crop` uses exactly the coordinates you pass.
It makes no assumption about what they mean, which allows it to accept any number of them.

:meth:`~sunpy.map.GenericMap.submap` treats its two arguments as opposite corners of a rectangle *in world coordinates*, then derives the other two corners before cropping.
It combines the supplied longitude and latitude values, so ``(-900", -900")`` and ``(0", -200")`` imply two additional corners: ``(-900", -200")`` and ``(0", -900")``.
It passes all four corners to :meth:`~ndcube.NDCube.crop`.

You can get the same behavior from :meth:`~ndcube.NDCube.crop` by passing all four corners yourself; :meth:`~sunpy.map.GenericMap.submap` fills in the two extra points for you.

For many maps, the pixel grid is aligned with the world coordinate axes: moving along a row changes longitude, and moving along a column changes latitude.
The rectangle defined by the two supplied corners is then aligned with the pixels. The two derived corners fall on its edges, so both methods return the same result.

On a rotated map
================

When the pixel grid is not aligned with the world coordinate axes, as in a rotated map, a rectangle bounded by straight lines in world coordinates is tilted relative to the pixel grid.
The smallest pixel rectangle containing it is larger than the tilted region itself.

The two supplied coordinates form one diagonal of the tilted rectangle; the derived corners form the other.
Each diagonal has its own pixel bounding box, and on a rotated map those boxes differ:

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

Drawn on the pixel grid, with the tilted rectangle in grey:

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

The grey rectangle is the requested region. ``bl`` and ``tr`` are the supplied coordinates, and the open markers show the corners derived by :meth:`~sunpy.map.GenericMap.submap`.

The crop is a narrow strip because its diagonal runs steeply across the pixel grid: it is tall but barely wide.
The submap also covers the other diagonal, which runs almost horizontally and extends the region sideways.
Both results have the same height because the first diagonal already spans the full extent in that direction.

Which one to use
================

Use :meth:`~sunpy.map.GenericMap.submap` when you want the region between two corners, interpreted as a rectangle on the sky.
Because the rectangle is defined in world coordinates, it covers the same part of the Sun regardless of the map's rotation, though the pixel cutout may be larger when the map is rotated.

Use :meth:`~ndcube.NDCube.crop` when you need the result to contain a set of coordinates, such as the corners of a tracked feature, without adding any implied corners.

To define a rectangle aligned with the pixel grid rather than the sky, pass pixel coordinates to :meth:`~sunpy.map.GenericMap.submap`.
It fills in the other corners on the pixel grid, so the result stays within the box defined by your two pixel positions.

In summary
==========

Given the same two coordinates, :meth:`~ndcube.NDCube.crop` uses a subset of the points used by :meth:`~sunpy.map.GenericMap.submap`.
The bounding box of a subset cannot be larger than that of the full set, so the region returned by :meth:`~ndcube.NDCube.crop` is always contained within the region returned by :meth:`~sunpy.map.GenericMap.submap`.

The methods return the same region when the pixel grid is aligned with the world coordinate axes, including after quarter turns.
Otherwise, :meth:`~ndcube.NDCube.crop` returns a smaller region. This includes rotations that are not multiples of 90 degrees and skewed grids with non-perpendicular axes.

This comparison applies only when both methods receive the same two coordinates.
If you give :meth:`~ndcube.NDCube.crop` more points, it returns the smallest region containing all of them, which may be larger.

See :ref:`sphx_glr_generated_gallery_map_submaps_and_cropping.py` and :ref:`sphx_glr_generated_gallery_map_cropping_a_map.py` for worked examples.
