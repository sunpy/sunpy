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
command only work once. If you were to call `~matplotlib.pyplot.show` again
after the above code is executed nothing happens. This confusing behavior
is something that the matplotlib devs get complaints about often and so this may change.
A discussion about this can be found `here
<http://stackoverflow.com/questions/5524858/matplotlib-show-doesnt-work-twice>`_.
Don't be confused by another command called `~matplotlib.pyplot.draw`.
This is only used while in interactive mode.

To turn on interactivity for pyplot use the command ::

    plt.ion()

In interactive mode, the plot will appear at the first `~matplotlib.pyplot.plot`
command and most commands will update the plot as you call them. Here is some
example code::

    plt.plot(range(10), range(10))
    plt.title("Simple Plot")

In this example, you'll see that the title appears right on the plot when you call it.
Note that in this case the `~matplotlib.pyplot.show` command is useless as the
plot shows up right when you create it. Also note that some commands will not
automatically update the plot and you have to use the `~matplotlib.pyplot.draw`
command. The following command ::

    plt.ioff()

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
 creates a figure then creates an axes
and populates the plot in ``ax``. With this method you
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

5. peek()
---------

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
plot window blocks and must be closed before doing anything else.

6. plot()
---------

For more advanced plotting the base SunPy objects also provide a `~sunpy.map.mapbase.GenericMap.plot` command.
This command is similar to the pyplot `~matplotlib.pyplot.plot` command in that
it will create a figure and axes object for you if you haven't already. It
returns a figure object and does not create a plot window. With the `~matplotlib.figure.Figure` object
in your hands you can reach in and grab the axes and therefore manipulate the plot.
Here is a simple example which outputs the same plot as we saw before. Click
on the link to see the code.

.. plot::

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

By default :ref:`map` uses the `wcsaxes <http://wcsaxes.readthedocs.org/>`_
package to improve the representation of world coordinates. In the
examples above the axes were normal matplotlib axes.
To create a custom `wcsaxes.WCSAxes` instance do the following ::

    fig = plt.figure()
    ax = plt.subplot(projection=smap.wcs)

when overplotting data and using wcsaxes you have to use the transform keyword
argument, also the native coordinate system of a `~wcsaxes.WCSAxes` is always
in degrees ::

    fig = plt.figure()
    ax = plt.subplot(projection=smap.wcs)

    smap.plot()
    ax.plot((100*u.arcsec).to(u.deg), (500*u.arcsec).to(u.deg),
            transform=ax.get_transform('world'))

Finally, here is a more complex example whose source code is available through
the link.

.. plot::

    from matplotlib import patches
    import astropy.units as u

    import sunpy.map
    import matplotlib.pyplot as plt
    import sunpy.data.sample

    smap = sunpy.map.Map(sunpy.data.sample.AIA_171_IMAGE)
    submap = smap.submap([-100-250, -100+250]*u.arcsec, [-400-250, -400+250]*u.arcsec)
    rect = patches.Rectangle([-100-250, -400-250], 500, 500, color = 'white', fill=False)

    fig = plt.figure()
    ax1 = fig.add_subplot(2,1,1)
    smap.plot()
    ax1.add_artist(rect)

    ax2 = fig.add_subplot(2,1,2)
    submap.plot()
    submap.draw_grid(grid_spacing=10*u.deg)
    ax2.set_title('submap')
    fig.subplots_adjust(hspace=0.4)

    plt.show()
