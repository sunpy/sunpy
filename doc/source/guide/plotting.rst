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
necessary to turn this feature on. This will be discussed below.
Pylab is another matplotlib usage scenario but it is essentially just pyplot with the
interactive capabilities turned on and numpy and matplotlib imported into the main
namespace.

2. Pyplot
---------
Here is an example of pyplot usage ::

    import matplotlib.pyplot as plt

    plt.plot(range(10), range(10))
    plt.title("A simple Plot")
    plt.show()

The show() command opens a plot on the screen and blocks execution until the plot window is closed. For
The show() command only work once. If you were to call show() on the plt object again
after the above code is executed nothing happens. Apparently, this confusing behavior
is something that the matplotlib devs get complaints about often and so this may change.
A discussion about this can be found `here
<http://stackoverflow.com/questions/5524858/matplotlib-show-doesnt-work-twice>`_.
Don't be confused by another command called draw(). This is only used while in interactive
mode.

To turn on interactivity for pyplot use the command ::

    plt.ion()

In interactive mode, the plot will appear at the first plot() command and most
commands will update the plot as you call them. Here is an example ::

    plt.plot(range(10), range(10))
    plt.title("Simple Plot")

In this example, you'll see that the title appears right on the plot when you call it.
Note that in this case the show command is useless as the plot shows up right when you
create it. Also note that some commands will not automatically update the plot and
you have to use the draw() command. The following command ::

    plt.ioff()

turns off interactivity.

3. Advanced Pyplot
------------------
If you need more fine-grained control over plots the recommended path is to use pyplot
and access the figures and axes objects. Here is an example ::

    import matplotlib.pyplot as plt
    import numpy as np

    x = np.arange(0, 10, 0.2)
    y = np.sin(x)

    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(x, y)
    ax.set_xlabel('x')

    plt.show()

Figure is the top-level container for all plot elements and axes is the top-level container
for a particular plot. So the above example, creates a figure then creates an axes
and then uses the pyplot plot() method to populate the plot in ax. With this method you
now have your hands on the ax
object so you can do things like change the labels on the x and y axes or add a legend.
In the previous section, pyplot took care of creating these
objects for you so you don't have to worry about creating them yourself.

4. SunPy Plotting Standards
---------------------------

To be consistent with matplotlib, SunPy has developed a standard plotting policy which
supports both simple and advanced matplotlib usage.

5. peek()
---------

For quick and easy access to a plot
all sunpy base objects (e.g. maps, spectra, lightcurves) define their own peek() command.
For example you can do the following ::

    import sunpy.map
    import sunpy.data.sample
    smap = sunpy.map.Map(sunpy.data.sample.EIT_195_IMAGE)
    smap.peek(draw_limb=True)

This creates a plot window with all axes defined, a plot title, and the image of the map
data defined by the contents of the map. As this command makes use of show(), in non-interactive
mode the plot window blocks and must be closed before doing anything else. This is meant as a
quick way to visualize the contents of a sunpy object you've created.

6. plot()
---------

For more advanced plotting the base sunpy objects also provide a plot() command. This
command is similar to the pyplot plot() command in that it will create a figure and axes
object for you if you haven't already. It returns a figure object and does not create a
plot window. With the figure object in your hands you can reach in and grab the axes
and therefore manipulate the plot. Here is an example of this at work ::

    import sunpy.map
    import sunpy.data.sample
    import matplotlib.pyplot as plt

    smap = sunpy.map.Map(sunpy.data.sample.EIT_195_IMAGE)
    smap.plot()
    smap.draw_limb()

    plt.show()

This output of this example is equivalent to one in the previous section. The `sunpy.map.Map.plot`
command is equivalent to the `~matplotlib.axes.Axes.imshow` command.
Similar to that command it will create a figure for you if you haven't created on yourself. For
advanced plotting you'll want to create it yourself. ::

    fig = plt.figure()
    ax = plt.subplot()

    smap.plot()
    plt.colorbar()
    ax.plot([-1000,1000], [0,0], color="white")

    plt.show()

The above will plot of line across the map. Using the fig.gca() command to get access to the
axes object most anything can be done to the plot and the plot can be displayed as usual
using the `~matplotlib.pyplot.show` command. Here is another example ::

    from matplotlib import patches
    fig = plt.figure()
    ax = plt.subplot()

    smap.plot()
    rect = patches.Rectangle([-350, -650], 500, 500, color = 'white', fill=False)
    ax.add_artist(rect)

    plt.show()

By default `~sunpy.map.Map` uses the `wcsaxes <http://wcsaxes.readthedocs.org/>`_
package to improve the representation of world coordinates on plots. In the
examples above the axes created is a normal matplotlib axes.
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

Finally, here is a more complex example::

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

The above example creates two side by side plots one with the overall view of the Sun
with a small area marked with a white box. That smaller view is then shown in the plot
below it. The spacing between the two plots is controlled by fig.subplots_adjust().

7. Plotting Keywords
--------------------

As mentioned before for Map `~matplotlib.pyplot.imshow()` does most of the heavy
lifting in the background while SunPy makes a number of choices for you so that
you don't have to. Changing these defaults
is made possible through some simple interfaces. Firstly you can pass any
`~matplotlib.pyplot.imshow()` keyword into
the plot command to override the defaults for that particular plot. The following example
changes the default AIA color table to use an inverse Grey color table::

    import sunpy.map
    import sunpy.data.sample
    import matplotlib.pyplot as plt
    smap = sunpy.map.Map(sunpy.data.sample.AIA_171_IMAGE)

    fig = plt.figure()
    ax = plt.subplot(1,1,1)
    smap.plot(cmap=plt.Greys_r)
    plt.show()

If you'd like to make this a permanent change you can access a number of settings under the
`plot_settings` property to make your changes for that map instance permanent.
 In the following example we change the normalization of the color table to a linear
one running from 5 to 100 (clipping everything above and below these values)::

    import sunpy.map
    import sunpy.data.sample
    import matplotlib.colors as colors
    smap = sunpy.map.Map(sunpy.data.sample.AIA_171_IMAGE)
    smap.plot_settings['norm'] = colors.Normalize(5, 1000)

    fig = plt.figure()
    ax = plt.subplot(1,1,1)
    smap.plot()
    plt.show()

8. Colormaps
------------

There are a number of color maps defined in SunPy which are used for data from
particular missions (e.g. SDO/AIA). The Map object chooses the appropriate colormap
on its own when you create it. The following example will show you all of the
colormaps available::

    import matplotlib.pyplot as plt
    import sunpy.cm

    # Access SunPy colormaps through matplotlib
    # You need to import sunpy.cm or sunpy.map for this to work.
    cmap = plt.get_cmap('sdoaia171')

    # Get a list of SunPy colormaps
    sunpy.cm.cmlist.keys()

    # you can also get a visual representation of all of the color tables
    sunpy.cm.show_colormaps()


.. image:: ../images/plotting_ex2.png

These can be used with the standard commands to change the colormap. So for
example if you wanted to plot an AIA image but use an EIT colormap, you would
do so as follows::

    import sunpy.map
    import sunpy.data.sample
    import matplotlib.pyplot as plt

    smap = sunpy.map.Map(sunpy.data.sample.AIA_171_IMAGE)
    cmap = plt.get_cmap('sohoeit171')

    fig = plt.figure()
    ax = plt.subplot(1,1,1)
    smap.plot(cmap=cmap)
    plt.show()
