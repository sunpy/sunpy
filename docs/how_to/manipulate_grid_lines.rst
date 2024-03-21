.. _how-to-manipulate-grid-lines-in-image-plots:

Manipulate grids lines when plotting Map with WCSAxes
=====================================================

Turning off the Helioprojective grid
------------------------------------

By default `sunpy.map.GenericMap.plot` draws a grid tracing the lines of helioprojective latitude and longitude, as below:

.. plot::
    :include-source:
    :context: close-figs

    import matplotlib.pyplot as plt
    import sunpy.map
    from sunpy.data.sample import AIA_171_IMAGE

    smap = sunpy.map.Map(AIA_171_IMAGE)

    fig = plt.figure()
    ax = fig.add_subplot(111, projection=smap)

    smap.plot(axes=ax)

    plt.show()

In some cases it may be desirable to remove this grid.
Since the underlying axes is a `~astropy.visualization.wcsaxes.WCSAxes`, you will need to access the ``grid`` method on the axes' object ``coord`` attribute:

.. plot::
    :include-source:
    :context: close-figs

    import matplotlib.pyplot as plt
    import sunpy.map
    from sunpy.data.sample import AIA_171_IMAGE

    smap = sunpy.map.Map(AIA_171_IMAGE)

    fig = plt.figure()
    ax = fig.add_subplot(111, projection=smap)

    smap.plot(axes=ax)
    ax.coords[0].grid(draw_grid=False)  # Disable grid for 1st axis
    ax.coords[1].grid(draw_grid=False)  # Disable grid for 2nd axis

    plt.show()

Changing the appearance of Heliographic grid
--------------------------------------------

`sunpy.map.GenericMap.draw_grid` allows users to overlay a Heliographic grid on their solar image plot.
This method does not explicitly provide many options for  the appearance of the grid lines.
Instead is accepts keyword arguments and transparently passes them onto the underlying infrastructure.
Therefore to change the width of the lines to say, 0.25, and their style to say, dotted, do the following:

.. plot::
    :include-source:
    :context: close-figs

    import matplotlib.pyplot as plt
    import sunpy.map
    from sunpy.data.sample import AIA_171_IMAGE

    smap = sunpy.map.Map(AIA_171_IMAGE)

    fig = plt.figure()
    ax = fig.add_subplot(111, projection=smap)

    smap.plot(axes=ax)
    smap.draw_grid(axes=ax, linewidth=0.25, linestyle="dotted")

    plt.show()
