.. _how-to-manipulate-grid-lines-in-image-plots:

Manipulate Grids Lines in Image Plots
=====================================

Turning Off Helioprojective Grid
--------------------------------

By default `sunpy.map.Map.plot` draws a grid tracing the lines of helioprojective
latitude and longitude, as below::

    >>> import matplotlib.pyplot as plt
    >>> import sunpy.map
    >>> from sunpy.data.sample import AIA_171_IMAGE

    >>> smap = sunpy.map.Map(AIA_171_IMAGE)

    >>> fig = plt.figure()
    >>> ax = fig.add_subplot(111, projection=smap)
    >>> smap.plot(ax)
    >>> plt.show()


In some cases it may be desirable to remove this grid. This can be done by
accesing the ``grid`` method on the axes' object ``coord`` attribute::

    >>> fig = plt.figure()
    >>> ax = fig.add_subplot(111, projection=smap)
    >>> smap.plot(ax)
    >>> ax.coords[0].grid(draw_grid=False)  # Disable grid for 1st axis
    >>> ax.coords[1].grid(draw_grid=False)  # Disable grid for 2nd axis
    >>> plt.show()

Changing Appearance of Heliocentric Grid
----------------------------------------

`sunpy.ma.Map.draw_grid` allows users to overlay a heliocentric grid on their
solar image plot::

    >>> fig = plt.figure()
    >>> ax = fig.add_subplot(111, projection=smap)
    >>> smap.plot(ax)

This method does not explicitly provide many options for manipulating
the appearance of the grid lines. Instead is accepts keyword arguments and
transparently passes them onto the underlying infrastructure.
Therefore to change the width of the lines to say, 0.25, and their style to
say, dotted, do the following::

    >>> fig = plt.figure()
    >>> ax = fig.add_subplot(111, projection=smap)
    >>> smap.plot(ax)
    >>> smap.draw_grid(ax, linewidth=0.25, linestyle="dotted")
    >>> plt.show()
