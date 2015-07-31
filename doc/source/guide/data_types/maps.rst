====
Maps
====

Maps in SunPy are are 2-dimensional data associated with a coordinate system.
In this guide, we will cover some of the basic functionality of maps.
Once you've read through this guide check out
the :doc:`/code_ref/map` for a more thorough look at SunPy maps.
There you can see what instruments are currently supported or you can access the
code reference for each instrument-specific map subclass.

Creating maps
-------------
To make things easy, SunPy can download several example files which are used
throughout the docs. These files have names like
`~sunpy.data.sample.AIA_171_IMAGE` and `~sunpy.data.sample.RHESSI_IMAGE`.
To create the sample `sunpy.map.sources.sdo.AIAMap` type the following into your
interactive Python shell::

    >>> import sunpy
    >>> import sunpy.map
    >>> import sunpy.data.sample
    >>> my_map = sunpy.map.Map(sunpy.data.sample.AIA_171_IMAGE)

If you have not downloaded the data already you should get an error and some
instruction on how to download the sample data.

The variable my_map is a :ref:`map` object. To create one from a
local FITS file try the following::

    >>> my_map = sunpy.map.Map('/mydirectory/mymap.fits')   # doctest: +SKIP

SunPy should automatically detects the type of file (e.g. FITS), what instrument it is
associated with (e.g. AIA, EIT, LASCO) and will automatically look in the
appropriate places for the FITS keywords it needs to interpret the coordinate
system. If the type of FITS file is not recognized then SunPy will try some
default FITS keywords and return a `~sunpy.map.GenericMap` but results
may vary. SunPy can also create maps from the jpg2000 files from
`helioviewer.org <http://helioviewer.org/>`_.

Creating Custom Maps
--------------------
It is also possible to create maps using custom data (e.g. from a simulation).
To do this you need to provide `~sunpy.map.map_factory.MapFactory` with both the data array as
well as some basic meta information. If no header is given then some default
values as assumed. Here is a simple example::

    >>> import numpy as np
    >>> data = np.arange(0,100).reshape(10,10)
    >>> header = {'cdelt1': 10, 'cdelt2': 10, 'telescop':'sunpy'}
    >>> my_map = sunpy.map.Map(data, header)

The keys in the header follows the `FITS standard <http://fits.gsfc.nasa.gov/fits_dictionary.html>`_.

Inspecting maps
---------------
A map contains a number of data-associated attributes. To get a quick look at
your map simply type::

    >>> my_map = sunpy.map.Map(sunpy.data.sample.AIA_171_IMAGE)
    >>> my_map   # doctest: +NORMALIZE_WHITESPACE
    SunPy AIAMap
    ---------
    Observatory:         SDO
    Instrument:  AIA 3
    Detector:    AIA
    Measurement:         171.0 Angstrom
    Wavelength:  171.0 Angstrom
    Obs Date:    2011-03-19 10:54:00
    dt:          1.999601 s
    Dimension:   [ 1024.  1024.] pix
    scale:               [ 2.4  2.4] arcsec / pix
    <BLANKLINE>
    array([[ 0.3125, -0.0625, -0.125 , ...,  0.625 , -0.625 ,  0.    ],
           [ 1.    ,  0.1875, -0.8125, ...,  0.625 , -0.625 ,  0.    ],
           [-1.1875,  0.375 , -0.5   , ..., -0.125 , -0.625 , -1.1875],
           ...,
           [-0.625 ,  0.0625, -0.3125, ...,  0.125 ,  0.125 ,  0.125 ],
           [ 0.5625,  0.0625,  0.5625, ..., -0.0625, -0.0625,  0.    ],
           [ 0.5   , -0.125 ,  0.4375, ...,  0.6875,  0.6875,  0.6875]])


This will show a representation of the data as well as some of its associated
attributes. A number of other attributes are also available, for example the
`~sunpy.map.GenericMap.date`, `~sunpy.map.GenericMap.exposure_time`,
`~sunpy.map.GenericMap.center`, `~sunpy.map.GenericMap.xrange`,
`~sunpy.map.GenericMap.yrange` and others (see `~sunpy.map.GenericMap`)::

    >>> map_date = my_map.date
    >>> map_exptime = my_map.exposure_time
    >>> map_center = my_map.center
    >>> map_xrange = my_map.xrange
    >>> map_yrange = my_map.yrange

To get a list of all of the attributes check the documentation by typing::

    >>> help(my_map)   # doctest: +SKIP

From SunPy version 0.6 many attributes return `~astropy.units.quantity.Quantity`
objects, please refer to the `Astropy documentation <http://astropy.readthedocs.org/en/latest/units/>`_.

The meta data for the map is accessed by ::

    >>> header = my_map.meta

This references the meta data dictionary with the header information as read
from the source file.

Getting at the data
-------------------
The data in a SunPy Map object is accessible through the
`~sunpy.map.GenericMap.data` attribute.  The data is implemented as a
NumPy `~numpy.ndarray`, so for example, to get
the 0th element in the array ::

    >>> my_map.data[0, 0]
    0.3125
    >>> my_map.data[0][0]
    0.3125

One important fact to remember is that the first
index is for the y direction while the second index is for the x direction.
For more information about indexing please refer to the
`Numpy documentation <http://www.scipy.org/Tentative_NumPy_Tutorial#head-864862d3f2bb4c32f04260fac61eb4ef34788c4c>`_.

Data attributes like `~numpy.ndarray.dtype` and
`~sunpy.map.GenericMap.dimensions` are accessible through
the SunPyGenericMap object ::

    >>> my_map.dimensions
    Pair(x=<Quantity 1024.0 pix>, y=<Quantity 1024.0 pix>)
    >>> my_map.dtype
    dtype('float64')

Here the dimensions attribute is similar to the `~numpy.ndarray.shape`
attribute, however returning an `~astropy.units.quantity.Quantity`.

If you'd like to use the data in a SunPy `~sunpy.map.GenericMap` object
elsewhere, you can use either of the following::

    >>> var = my_map.data
    >>> var = my_map.data.copy()

Python makes use of pointers so if you want to alter the data and keep the
original data in the map intact make sure to copy it.

Some basic statistical functions on the data array are also passed through to Map
objects::

    >>> my_map.min()
    -2.0
    >>> my_map.max()
    9429.125
    >>> my_map.mean()
    235.91531443595886

but you can also access all the other `~numpy.ndarray` functions and attributes
but accessing the data array directly. For example::

    >>> my_map.data.std()
    292.43424704677756

Plotting
--------
As is true of all of the SunPy data objects, the SunPy `~sunpy.map.GenericMap`
object (and all of its instrument-specific sub-classes) has its
own built-in plot methods so that it is easy to quickly view your map.
To create a plot just type::

    >>> my_map.peek()   # doctest: +SKIP

This will open a matplotlib plot on your screen.
In addition, to enable users to modify the plot it is possible to grab the
matplotlib axes object by using the `~sunpy.map.GenericMap.plot()` command.
This makes it possible to use the SunPy plot as the foundation for a
more complicated figure. For a bit more information about this and some
examples see :ref:`plotting`.

.. note::

   If the `wcsaxes <http://wcsaxes.readthedocs.org/en/latest/>`_ package is not
   installed the `~sunpy.map.GenericMap.plot()` and `~sunpy.map.GenericMap.peek()`
   methods assume that the data is not rotated,
   i.e. the solar y axis is oriented with the columns of the array. If this condition
   is not met, when the map is plotted a warning will be issued. You can create
   an oriented map by using `~sunpy.map.GenericMap.rotate()` before you plot the Map.

Plotting Keywords
*****************

For Map `~matplotlib.pyplot.imshow` does most of the heavy
lifting in the background while SunPy makes a number of choices for you so that
you don't have to (e.g. colortable, plot title). Changing these defaults
is made possible through two simple interfaces. You can pass any
`~matplotlib.pyplot.imshow` keyword into
the plot command to override the defaults for that particular plot. The following
plot changes the default AIA color table to use an inverse Grey color table.

.. plot::
    :include-source:

    import sunpy.map
    import sunpy.data.sample
    import matplotlib.pyplot as plt
    smap = sunpy.map.Map(sunpy.data.sample.AIA_171_IMAGE)
    fig = plt.figure()
    smap.plot(cmap=plt.cm.Greys_r)
    plt.colorbar()
    plt.show()

You can view or make changes to the default settings through the `~sunpy.map.GenericMap.plot_settings`
dictionary. In the following example we change the title of the plot by changing the
`~sunpy.map.GenericMap.plot_settings` property.

.. plot::
    :include-source:

    import sunpy.map
    import sunpy.data.sample
    import matplotlib.pyplot as plt
    smap = sunpy.map.Map(sunpy.data.sample.AIA_171_IMAGE)
    smap.plot_settings['title'] = "My Second SunPy Plot"
    smap.plot_settings['cmap'] = plt.cm.Blues_r
    fig = plt.figure()
    smap.plot()
    plt.colorbar()
    plt.show()


Colormaps and Normalization
***************************

Image data is generally shown in false color in order to better identify it or
to better visualize structures in the image. Matplotlib handles this colormapping
process through the `~matplotlib.colors` module. This process involves two steps:
the data array is first mapped onto the range 0-1 using an instance of
`~matplotlib.colors.Normalize` or a subclass; then this number is mapped to a
color using an instance of a subclass of a `~matplotlib.colors.colormap`.

SunPy provides the colormaps for each mission as defined by the mission teams.
The Map object chooses the appropriate colormap for you when it is created as
long as it recognizes the instrument. To see what colormaps are available::

    >>> import sunpy.cm
    >>> sunpy.cm.cmlist.keys()   # doctest: +NORMALIZE_WHITESPACE
    ['sohoeit304', 'sdoaia211', 'sohoeit195', 'trace1600', 'sdoaia94',
     'trace284', 'trace1216', 'sdoaia304', 'trace1700', 'yohkohsxtal',
     'trace195', 'sdoaia335', 'sdoaia1600', 'traceWL', 'stereohi2',
     'sdoaia193', 'stereohi1', 'rhessi', 'trace171', 'trace1550',
     'sohoeit284', 'stereocor2', 'hmimag', 'stereocor1', 'sdoaia1700',
     'yohkohsxtwh', 'sohoeit171', 'hinodexrt', 'sdoaia131', 'sdoaia171',
     'hinodesotintensity', 'sdoaia4500', 'soholasco3', 'soholasco2']

The SunPy colormaps are registered with matplotlib so you can grab them like
you would any other colormap::

    >>> import matplotlib.pyplot as plt   # doctest: +SKIP
    >>> import sunpy.cm

You need to import sunpy.cm or sunpy.map for this to work::

    >>> cmap = plt.get_cmap('sdoaia171')   # doctest: +SKIP


The following plot shows off all of the colormaps.

.. plot::

    import matplotlib.pyplot as plt
    import sunpy.cm
    sunpy.cm.show_colormaps()

These can be used with the standard commands to change the colormap. So for
example if you wanted to plot an AIA image but use an EIT colormap, you would
do so as follows.

.. plot::
    :include-source:

    import sunpy.map
    import sunpy.data.sample
    import matplotlib.pyplot as plt

    smap = sunpy.map.Map(sunpy.data.sample.AIA_171_IMAGE)
    cmap = plt.get_cmap('sohoeit171')

    fig = plt.figure()
    ax = plt.subplot(1,1,1)
    smap.plot(cmap=cmap)
    plt.colorbar()
    plt.show()

or you can just change the colormap for the map itself as follows::

    >>> smap.plot_settings['cmap'] = plt.get_cmap('sohoeit171')   # doctest: +SKIP

The normalization is also set automatically and is chosen so that all the
data from minimum to maximum is displayed as best as possible for most cases.
This means that it is never necessary to touch the data such as applying a function
such sqrt or log to the data to make your plot look good.
There are many normalizations available from matplotlib such as `~matplotlib.colors.Lognorm`. Other
`more exotic normalizations <http://docs.astropy.org/en/stable/visualization/index.html>`_ are also
made available from Astropy.  Just like the colormap the default normalization
can be changed through the plot_settings dictionary or directly for the individual
plot by passing a keyword argument. The following example shows the difference between
a linear and logarithmic normalization on an AIA image.

.. plot::
    :include-source:

    import sunpy.map
    import sunpy.data.sample
    import matplotlib.pyplot as plt
    import matplotlib.colors as colors

    smap = sunpy.map.Map(sunpy.data.sample.AIA_171_IMAGE)

    fig = plt.figure()
    ax1 = fig.add_subplot(2,1,1)
    smap.plot(norm=colors.Normalize())
    plt.colorbar()
    ax2 = fig.add_subplot(2,1,2)
    smap.plot(norm=colors.LogNorm())
    fig.subplots_adjust(hspace=0.4)
    plt.colorbar()
    plt.show()

Note how the color in the colorbar does not change since these two maps share
the same colormap while the data values associated with each color do because
the normalization is different.

Masking and Clipping Data
-------------------------
It is often necessary for the purposes of display or otherwise to ignore certain
data in an image. For example large data value could be due to
cosmic ray hits and should be ignored. The most straightforward way to ignore
this kind of data in plots without altering the data is to clip it. This can be achieved
very easily when initializing the normalization variable. For example::

    >>> import matplotlib.colors as colors
    >>> norm = colors.Normalize(vmin=smap.min(), vmax=smap.mean() + 3 *smap.std())   # doctest: +SKIP

This clips out many of the brightest pixels. If you'd like to see what areas of
your images got clipped set the following values::

    >>> cmap = cmap.plot_settings['cmap']   # doctest: +SKIP
    >>> cmap.set_over('red', 1.0)   # doctest: +SKIP
    >>> cmap.set_under('green', 1.0)   # doctest: +SKIP

This will color the areas above and below in red and green respectively
(similar to this `example <http://matplotlib.org/examples/pylab_examples/image_masked.html>`_).
You can use the following colorbar command to display these choices::

    >>> plt.colorbar(extend='both')   # doctest: +SKIP

Here is an example of this put to use on an AIA image. If you see how the image
displays by default you'll see that it does not look that pretty. This is
because the image contains some negative values which are throwing off the
normalization.

.. plot::

    import sunpy.map
    import matplotlib.pyplot as plt
    import sunpy.data.sample
    smap = sunpy.map.Map(sunpy.data.sample.AIA_94_CUTOUT)
    txt = "min={min}, max={max}, $\mu$={mean}, $\sigma$={std}".format(min=int(smap.min()),
                                                                      max=int(smap.max()),
                                                                      mean=int(smap.mean()),
                                                                      std=int(smap.std()))
    plt.text(-1100, 0, txt, color='white')
    smap.plot()
    plt.colorbar()
    plt.show()

In order to fix this we need to adjust our normalization to not display negative
values. We can also brighten the image by clipping the high values though this
will mean that the bright regions look 'saturated'. This is achieved in the following plot.

.. plot::
    :include-source:

    import sunpy.map
    import matplotlib.pyplot as plt
    import matplotlib.colors as colors
    import sunpy.data.sample
    smap = sunpy.map.Map(sunpy.data.sample.AIA_94_CUTOUT)
    cmap = smap.plot_settings['cmap']
    cmap.set_over('blue', 1.0)
    cmap.set_under('purple', 1.0)
    norm = colors.Normalize(vmin=0, vmax=smap.mean() + 5 * smap.std())
    smap.plot(norm=norm)
    plt.colorbar(extend='both')
    plt.show()

Another method to ignore bad data is to mask the data. A mask is a boolean
array and so can give you much more fine-grained control over what is not being
displayed.  A `~numpy.ma.MaskedArray`
is a subclass of a numpy array so it has all of the same properties with the
addition of an associated boolean array which holds the mask.

.. the following is a good example which could be fixed and added later
.. The following plot achieves the same goal as above but using a mask instead of clipping.

..    import sunpy.map
    import matplotlib.pyplot as plt
    import matplotlib.colors as colors
    cmap = smap.plot_settings['cmap']
    cmap.set_bad('blue', 1.0)
    smap = sunpy.map.Map('/Users/schriste/Downloads/old downloads/foxsi_ar_data/ssw_cutout_20121030_153001_AIA_94_.fts')
    smap.mask =
    smap.plot()
    plt.colorbar(extend='both')
    plt.show()

.. Hinode XRT image. By inspecting the maximum versus the mean and standard deviation, it is clear that there are some overly bright pixels. This is likely due to cosmic ray hits which is throwing off the default plot making it too dark to see the solar emission.

.. .. plot::

..    import sunpy.map
    import matplotlib.pyplot as plt
    smap = sunpy.map.Map('/Users/schriste/Desktop/sunpy_test_img/XRT20141211_184221.9.fits')
    fig = plt.figure()
    smap.plot()
    txt = "min={min}, max={max}, $\mu$={mean}, $\sigma$={std}".format(min=int(smap.min()),
                                                                      max=int(smap.max()),
                                                                      mean=int(smap.mean()),
                                                                      std=int(smap.std()))
    plt.text(-600, 1500, txt, color='white')
    plt.colorbar()
    plt.show()

.. Let's address this by clipping the largest values (in this case everything above 3 sigma). The following plot shows the result of this operation.

.. .. plot::

..     import sunpy.map
    import matplotlib.pyplot as plt
    import matplotlib.colors as colors
    cmap = smap.plot_settings['cmap']
    cmap.set_over('green', 1.0)
    cmap.set_under('purple', 1.0)
    norm = colors.Normalize(vmin=smap.min(), vmax=smap.mean() + 3 *smap.std())
    smap = sunpy.map.Map('/Users/schriste/Desktop/sunpy_test_img/XRT20141211_184221.9.fits')
    smap.plot(norm=norm)
    plt.colorbar(extend='both')
    plt.show()

.. This makes it very visible that there are a number of hot pixels mostly concentrated in the upper half of this image. Now let's address this problem with masking instead of clipping.

.. .. plot::

..     import sunpy.map
    import matplotlib.pyplot as plt
    import matplotlib.colors as colors
    import numpy.ma
    smap = sunpy.map.Map('/Users/schriste/Desktop/sunpy_test_img/XRT20141211_184221.9.fits')
    cmap = smap.plot_settings['cmap']
    cmap.set_bad('blue', 1.0)
    smap.data = numpy.ma.masked_greater(smap.data, smap.mean() + 3 *smap.std())
    txt = "min={min}, max={max}, $\mu$={mean}, $\sigma$={std}".format(min=int(smap.min()),
                                                                      max=int(smap.max()),
                                                                      mean=int(smap.mean()),
                                                                      std=int(smap.std()))
    plt.text(-600, 1500, txt, color='white')
    norm = colors.Normalize()
    smap.plot(norm = norm)
    plt.colorbar(extend='both')

.. This plot shows a very similar effect to clipping but note that the array properties such as max and min have changed. That's because numpy is now ignoring those masked values. With a masked array
.. (compared to clipping) we can go ahead and make more detailed masking operations so that we are not masking the emission from the bright solar sources. The next plot masks only those bright pixels in the upper area of the plot leaving the bright solar sources which are concentrated in the lower part of the plot intact.

.. .. plot::

..     import sunpy.map
    import matplotlib.pyplot as plt
    import matplotlib.colors as colors
    import numpy.ma
    file = '/Users/schriste/Downloads/old downloads/foxsi_ar_data/sXRT20141211_184221.9.fits'
    smap = sunpy.map.Map(file)
    cmap = smap.plot_settings['cmap']
    cmap.set_bad('blue', 1.0)
    smap.data = numpy.ma.masked_greater(smap.data, smap.mean() + 3 *smap.std())
    smap.data.mask[0:250,:] = False
    txt = "min={min}, max={max}, $\mu$={mean}, $\sigma$={std}".format(min=int(smap.min()),
                                                                      max=int(smap.max()),
                                                                      mean=int(smap.mean()),
                                                                      std=int(smap.std()))
    plt.text(-600, 1500, txt, color='white')
    norm = colors.Normalize()
    smap.plot(norm = norm)
    plt.colorbar(extend='both')


Composite Maps and Overlaying Maps
----------------------------------

The `Map()` method described above can also handle a list of maps. If a series of maps
are supplied as inputs, `Map()` will return a list of maps as the output.  However,
if the 'composite' keyword is set to True, then a `~sunpy.map.CompositeMap` object is
returned.  This is useful if the maps are of a different type (e.g. different
instruments).  For example, to create a simple composite map::

    >>> my_maps = sunpy.map.Map(sunpy.data.sample.EIT_195_IMAGE, sunpy.data.sample.RHESSI_IMAGE, composite=True)

A `~sunpy.map.CompositeMap` is different from a regular SunPy `~sunpy.map.GenericMap` object and therefore
different associated methods. To list which maps are part of your composite map use::

    >>> my_maps.list_maps()
    [<class 'sunpy.map.sources.soho.EITMap'>, <class 'sunpy.map.sources.rhessi.RHESSIMap'>]

The following code adds a new map (which must be instantiated first), sets
its transparency to 25%, turns on contours from 50% to 90% for the second
map, and then plots the result.

.. plot::
    :include-source:

    import sunpy.data.sample
    import sunpy.map
    import matplotlib.pyplot as plt
    my_maps = sunpy.map.Map(sunpy.data.sample.EIT_195_IMAGE, sunpy.data.sample.RHESSI_IMAGE, composite=True)
    my_maps.add_map(sunpy.map.Map(sunpy.data.sample.AIA_171_IMAGE))
    my_maps.set_alpha(2, 0.5)
    my_maps.set_levels(1, [50, 60, 70, 80, 90], percent = True)
    my_maps.plot()
    plt.show()

This is not a particularly pretty plot but it shows what SunPy can do!

Working with your map
---------------------
Part of the philosophy of the map object is to provide most of the basic
functionality that a scientist would want therefore a map also contains a number
of map-specific methods such as resizing a map or grabbing a subview. To get
a list of the methods available for a map type::

    >>> help(my_map)   # doctest: +SKIP

and check out the methods section!

Mapcubes
--------
A `~sunpy.map.MapCube` is an ordered list of maps.  By default, the maps are ordered by
their observation date, from earlier maps to later maps. A `~sunpy.map.MapCube` can be
created by supplying multiple existing maps::

    >>> map1 = sunpy.map.Map(sunpy.data.sample.AIA_171_IMAGE)
    >>> map2 = sunpy.map.Map(sunpy.data.sample.EIT_195_IMAGE)
    >>> mc = sunpy.map.Map([map1, map2], cube=True)

or by providing a directory full of image files::

    >>> mc = sunpy.map.Map('path/to/my/files/*.fits', cube=True)   #  doctest: +SKIP

The earliest map in the mapcube can be accessed by simply indexing the maps
list::

    >>> mc.maps[0]   # doctest: +SKIP

Mapcubes can hold maps that have different shapes.  To test if all the
maps in a `~sunpy.map.MapCube` have the same shape::

    >>> mc.all_maps_same_shape()
    True

It is often useful to return the image data in a `~sunpy.map.MapCube` as a single
three dimensional Numpy `~numpy.ndarray`::

    >>> mc.as_array()   # doctest: +SKIP

Note that an array is returned only if all the maps have the same
shape.  If this is not true, an error (ValueError) is returned.  If all the
maps have nx pixels in the x-direction, and ny pixels in the y-direction,
and there are n maps in the mapcube, the `~numpy.ndarray` array that is
returned has shape (ny, nx, n).  The data of the first map in the `~sunpy.map.MapCube`
appears in the `~numpy.ndarray` in position ``[:, :, 0]``, the data of second map in
position ``[:, :, 1]``, and so on.  The order of maps in the `~sunpy.map.MapCube` is
reproduced in the returned `~numpy.ndarray`.

The meta data from each map can be obtained using::

    >>> mc.all_meta()   # doctest: +SKIP

This returns a list of map meta objects that have the same order as
the maps in the `~sunpy.map.MapCube`.

Coalignment of Mapcubes
-----------------------
A typical data preparation step when dealing with time series of images is to
coalign images taken at different times so that features in different images
remain in the same place.  A common approach to this problem is
to take a representative template that contains the features you are interested
in, and match that to your images.  The location of the best match tells you
where the template is in your image.  The images are then shifted to the
location of the best match.  This aligns your images to the position of the
features in your representative template.

SunPy provides a function to coalign the maps inside the `~sunpy.map.MapCube`.
The implementation of this functionality requires the installation of the
scikit-image library, a commonly used image processing library.
To coalign a `~sunpy.map.MapCube`, simply import
the function and apply it to your `~sunpy.map.MapCube`::

    >>> from sunpy.image.coalignment import mapcube_coalign_by_match_template
    >>> coaligned = mapcube_coalign_by_match_template(mc)

This will return a new `~sunpy.map.MapCube`, coaligned to a template extracted from the
center of the first map in the `~sunpy.map.MapCube`, with the map dimensions clipped as
required.  The coalignment algorithm provides many more options for handling
the coalignment of `~sunpy.map.MapCube` type::

    >>> help(mapcube_coalign_by_match_template)   # doctest: +SKIP

for a full list of options and functionality.

If you just want to calculate the shifts required to compensate for solar
rotation relative to the first map in the `~sunpy.map.MapCube` without applying them, use::

    >>> from sunpy.image.coalignment import calculate_match_template_shift
    >>> shifts = calculate_match_template_shift(mc)

This is the function used to calculate the shifts in `~sunpy.map.MapCube` coalignment
function above.  Please see `~sunpy.image.coalignment.calculate_match_template_shift` to learn more about its features.
Shifts calculated using calculate_match_template_shift can be passed directly
to the coalignment function.


Compensating for solar rotation in Mapcubes
-------------------------------------------
Often a set of solar image data consists of fixing the pointing of a
field of view for some time and observing.  Features on the Sun will
rotate according to the Sun's rotation.

A typical data preparation step when dealing with time series of these
types of images is to shift the images so that features do not appear
to move across the field of view.  This requires taking in to account
the rotation of the Sun.  The Sun rotates differentially, depending on
latitude, with features at the equator moving faster than features at
the poles.

SunPy provides a function to shift images in `~sunpy.map.MapCube` following solar
rotation.  This function shifts an image according to the solar
differential rotation calculated at the latitude of the center of the
field of view.  The image is not *differentially* rotated.  This
function is useful for de-rotating images when the effects of
differential rotation in the `~sunpy.map.MapCube` can be ignored (for example, if
the spatial extent of the image is small, or when the duration of the
`~sunpy.map.MapCube` is small; deciding on what 'small' means depends on your
application).

To apply this form of solar derotation to a `~sunpy.map.MapCube`, simply import the
function and apply it to your `~sunpy.map.MapCube`::

    >>> from sunpy.physics.transforms.solar_rotation import mapcube_solar_derotate
    >>> derotated = mapcube_solar_derotate(mc)

For more info see `~sunpy.physics.transforms.solar_rotation.mapcube_solar_derotate`.

If you just want to calculate the shifts required to compensate for solar
rotation relative to the first map in the `~sunpy.map.MapCube` without applying them, use::

    >>> from sunpy.physics.transforms.solar_rotation import calculate_solar_rotate_shift
    >>> shifts = calculate_solar_rotate_shift(mc)

Please consult the docstring of the `~sunpy.image.coalignment.mapcube_coalign_by_match_template` function in order to learn about
the features of this function.
