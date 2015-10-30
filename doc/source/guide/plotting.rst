.. _plotting:

-----------------
Plotting in SunPy
-----------------

SunPy makes use of `matplotlib <http://matplotlib.org/>`_ for all of its plotting needs
as such it tries to follow the matplotlib plotting philosophy.
It is therefore useful to go over how matplotlib works as background.

1. Matplotlib Tutorial
----------------------
The tutorial provided here is a summary of one that can be found in the `matplotlib
usage documentation <http://matplotlib.org/faq/usage_faq.html>`_.

Matplotlib provides two main pathways for plotting. One is meant for interactive use
(e.g. command-line) and the other for non-interactive use (e.g. modules). It is important
to recognize though that the interactive-use pathway (referred to as pyplot) just
provides shortcuts for doing many of the more advanced non-interactive functions in the
background. It is therefore possible to switch between the two as necessary and
it is possible to use pyplot in a non-interactive way. In this manner pyplot
is just a shortcut to making it quicker to set up plot axes and figures.
In order to get access to the full interactive capabilities of pyplot it is
necessary to turn this feature on.
Pylab is another matplotlib usage scenario but it is essentially just pyplot with the
interactive capabilities turned on and numpy and matplotlib imported into the main
namespace.

2. Pyplot
---------
Here is a simple example of pyplot usage.

.. plot::
    :include-source:

    import matplotlib.pyplot as plt
    plt.plot(range(10), range(10))
    plt.title("A simple Plot")
    plt.show()

The `~matplotlib.pyplot.show` command opens a plot on the screen and blocks
execution until the plot window is closed. The `~matplotlib.pyplot.show`
command only works once. If you were to call `~matplotlib.pyplot.show` again
after the above code is executed nothing happens. This confusing behavior
is something that the matplotlib devs get complaints about often and so this may change.
A discussion about this can be found `here
<http://stackoverflow.com/questions/5524858/matplotlib-show-doesnt-work-twice>`_.
Don't be confused by another command called `~matplotlib.pyplot.draw`.
This is only used while in interactive mode.

To turn on interactivity for pyplot use the command ::

    >>> plt.ion()   # doctest: +SKIP

In interactive mode, the plot will appear at the first `~matplotlib.pyplot.plot`
command and most commands will update the plot as you call them. Here is some
example code::

    >>> plt.plot(range(10), range(10))   # doctest: +SKIP
    >>> plt.title("Simple Plot")   # doctest: +SKIP

In this example, you'll see that the title appears right on the plot when you call it.
Note that in this case the `~matplotlib.pyplot.show` command is useless as the
plot shows up right when you create it. Also note that some commands will not
automatically update the plot and you have to use the `~matplotlib.pyplot.draw`
command. The following command ::

    >>> plt.ioff()   # doctest: +SKIP

turns off interactivity.

3. Advanced Pyplot
------------------
If you need more fine-grained control over plots the recommended path is to use pyplot
and access the figures and axes objects. This is shown in the following example.

.. plot::
    :include-source:

    import matplotlib.pyplot as plt
    import numpy as np
    x = np.arange(0, 10, 0.2)
    y = np.sin(x)
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(x, y)
    ax.set_xlabel('x')
    plt.show()

In matplotlib, `~matplotlib.figure.Figure` is the top-level container for all plot elements and
`~matplotlib.axes.Axes` is the top-level container for a particular plot. So the above example,
creates a figure then creates an axes and populates the plot in ``ax``. With this method you
now have your hands on the `~matplotlib.axes.Axes` object so you can do things
like change the labels on the x and y axes or add a legend.
In the previous section, pyplot took care of creating these
objects for you so you didn't have to worry about creating them yourself.

4. SunPy Plotting Standards
---------------------------

To be consistent with matplotlib, SunPy has developed a standard plotting policy
which supports both simple and advanced matplotlib usage. The following examples
focus on the map object but they should be applicable across all of the data
objects.

4.1 peek()
----------

For quick and easy access to a plot
all SunPy base objects (e.g. maps, spectra, lightcurves) define their own
`~sunpy.map.mapbase.GenericMap.peek` command which will create a plot for you and show it without you having to deal
with any matplotlib setup. This is so that it is easy to take a quick look at
your data. For example you can make the following plot.

.. plot::
    :include-source:

    import sunpy.map
    import sunpy.data.sample
    smap = sunpy.map.Map(sunpy.data.sample.AIA_171_IMAGE)
    smap.peek(draw_limb=True)

This creates a plot window with all axes defined, a plot title, and the image of
the map data defined by the contents of the map. In non-interactive mode the
plot window blocks the command line terminal and must be closed before doing anything else.

4.2 plot()
----------

For more advanced plotting the base SunPy objects also provide a `~sunpy.map.mapbase.GenericMap.plot` command.
This command is similar to the pyplot `~matplotlib.pyplot.plot` command in that
it will create a figure and axes object for you if you haven't already. It
returns a figure object and does not create a plot window. With the `~matplotlib.figure.Figure` object
in your hands you can reach in and grab the axes and therefore manipulate the plot.
Here is a simple example which outputs the same plot as we saw before:

.. plot::
    :include-source:

    import sunpy.map
    import sunpy.data.sample
    import matplotlib.pyplot as plt
    smap = sunpy.map.Map(sunpy.data.sample.AIA_171_IMAGE)
    smap.plot()
    smap.draw_limb()
    plt.colorbar()
    plt.show()

For more advanced plotting you'll want to create the `~matplotlib.figure.Figure` object yourself.
The following example plot shows how to add a rectangle to a plot to, for example,
highlight a region of interest, and change the plot title.

.. plot::
    :include-source:

    import sunpy.map
    import sunpy.data.sample
    import matplotlib.pyplot as plt
    from matplotlib import patches
    smap = sunpy.map.Map(sunpy.data.sample.AIA_171_IMAGE)

    fig = plt.figure()
    ax = plt.subplot()

    smap.plot()
    rect = patches.Rectangle([-350, -650], 500, 500, color = 'white', fill=False)
    ax.set_title('My customized plot')
    ax.add_artist(rect)
    plt.colorbar()
    plt.show()


Plotting Maps with wcsaxes
--------------------------

By default :ref:map checks if the `wcsaxes <http://wcsaxes.readthedocs.org/>`_ 
package has been installed. If it is installed, 
then `wcsaxes` is used to improve the representation of world coordinates,
and calling ~sunpy.map.GenericMap.plot or~sunpy.map.GenericMap.peek() will use 
wcsaxes for plotting. Unless a standard `matplotlib.axes.Axes` object is created.

To explicitly create a `wcsaxes.WCSAxes` instance do the following ::

    >>> fig = plt.figure()   # doctest: +SKIP
    >>> ax = plt.subplot(projection=smap.wcs)   # doctest: +SKIP

when plotting on a `~wcsaxes.WCSAxes` axes, it will by default plot in pixel 
coordinates, you can override this behavior and plot in 'world' coordinates
by getting the transformation from the axes with ``ax.get_transform('world')``.
Note: World coordinates are always in **degrees** so you will have to convert 
to degrees.::

    >>> smap.plot()   # doctest: +SKIP
    >>> ax.plot((100*u.arcsec).to(u.deg), (500*u.arcsec).to(u.deg),
    ...         transform=ax.get_transform('world'))   # doctest: +SKIP

Finally, here is a more complex example using SunPy maps, wcsaxes and Astropy 
units to plot a AIA image and a zoomed in view of an active region.

.. plot::
    :include-source:

    import matplotlib.pyplot as plt
    from matplotlib import patches
    import astropy.units as u

    import sunpy.map
    import sunpy.data.sample


    # Define a region of interest
    l = 250*u.arcsec
    x0 = -100*u.arcsec
    y0 = -400*u.arcsec

    # Create a SunPy Map, and a second submap over the region of interest.
    smap = sunpy.map.Map(sunpy.data.sample.AIA_171_IMAGE)
    submap = smap.submap(u.Quantity([x0-l, x0+l]), u.Quantity([y0-l, y0+l]))



    # Create a new matplotlib figure, larger than default.
    fig = plt.figure(figsize=(5,12))

    # Add a first Axis, using the WCS from the map.
    ax1 = fig.add_subplot(2,1,1, projection=smap.wcs)

    # Plot the Map on the axes with default settings.
    smap.plot()

    # Define a region to highlight with a box
    # We have to convert the region of interest to degress, and then get the raw values.
    bottom_left = u.Quantity([x0-l, y0-l]).to(u.deg).value
    l2 = (l*2).to(u.deg).value

    # create the rectangle, we use the world transformation to plot in physical units.
    rect = patches.Rectangle(bottom_left, l2, l2, color='white', fill=False,
                             transform=ax1.get_transform('world'))
                         
    # Add the rectangle to the plot.
    ax1.add_artist(rect)



    # Create a second axis on the plot.
    ax2 = fig.add_subplot(2,1,2, projection=submap.wcs)

    submap.plot()

    # Add a overlay grid.
    submap.draw_grid(grid_spacing=10*u.deg)

    # Change the title.
    ax2.set_title('Zoomed View')


    plt.show()
