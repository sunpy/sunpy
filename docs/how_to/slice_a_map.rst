.. _sunpy-how-to-slice-a-map:

*****************************
Slice a `~sunpy.map.Map`
*****************************

A `~sunpy.map.Map` can be sliced with array indices.
The result is a new map of the same type, whose metadata describes the region that was sliced out, so it plots and transforms correctly on its own.

Indices are in array order, ``map[rows, columns]``, which is the opposite order to the ``(x, y)`` ordering used by world coordinates.

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

Any mask and uncertainty are sliced alongside the data.

Slicing cannot reduce the number of dimensions of a map, so an integer index raises an error.
Use a length-one slice instead.

To cut out a region given in world coordinates rather than array indices, use :meth:`~sunpy.map.GenericMap.submap` or :meth:`ndcube.NDCube.crop` instead.
See :ref:`sphx_glr_generated_gallery_map_submaps_and_cropping.py` and :ref:`sphx_glr_generated_gallery_map_cropping_a_map.py` for examples, and :ref:`sunpy-topic-guide-submap-vs-crop` for a comparison of the two methods.
