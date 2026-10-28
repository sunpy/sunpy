.. _sunpy-how-to-slice-a-map:

*****************************
Slice a `~sunpy.map.Map`
*****************************

A `~sunpy.map.Map` can be sliced using array indices.
Slicing returns a new map of the same type, with metadata that describes the selected region so it plots and transforms correctly on its own.

Indices follow array order, ``map[rows, columns]``. World coordinates use the opposite order, ``(x, y)``.

.. plot::
    :include-source:
    :context: close-figs

    import matplotlib.pyplot as plt
    from matplotlib.patches import Rectangle

    import astropy.units as u

    import sunpy.data.sample
    import sunpy.map

    aia_map = sunpy.map.Map(sunpy.data.sample.AIA_304_IMAGE)
    cutout = aia_map[260:460, 680:880]

    fig = plt.figure(figsize=(10, 5))

    ax1 = fig.add_subplot(121, projection=aia_map)
    aia_map.plot(axes=ax1, clip_interval=(1, 99.99)*u.percent)
    ax1.add_patch(Rectangle((680, 260), 200, 200, edgecolor="cyan", facecolor="none"))

    ax2 = fig.add_subplot(122, projection=cutout)
    cutout.plot(axes=ax2, clip_interval=(1, 99.99)*u.percent)

    plt.show()

The slice includes any mask and uncertainty associated with the data.

Slicing cannot reduce a map's number of dimensions, so using an integer index raises an error.
Use a slice of length one instead.

To cut out a region using world coordinates, use :meth:`~sunpy.map.GenericMap.submap` or :meth:`ndcube.NDCube.crop`.
See :ref:`sphx_glr_generated_gallery_map_creating_a_submap.py` and :ref:`sphx_glr_generated_gallery_map_cropping_a_map.py` for examples.
The :ref:`sunpy-topic-guide-submap-vs-crop` topic guide explains how these two methods differ.
