****
Maps
****

Maps in SunPy are are 2-dimensional data associated with a coordinate system. In
this guide, we will cover some of the basic functionality of maps. Once you've
read through this guide check out :doc:`/code_ref/map` for a more thorough look
at SunPy maps. There you can see what instruments are currently supported or you
can access the code reference for each instrument-specific map subclass.

Creating maps
=============
To make things easy, SunPy can download several example files which are used
throughout the docs. These files have names like
``sunpy.data.sample.AIA_171_IMAGE`` and ``sunpy.data.sample.RHESSI_IMAGE``. To
create a `~sunpy.map.Map` from the the sample AIA image
type the following into your Python shell::

    >>> import sunpy
    >>> import sunpy.map
    >>> import sunpy.data.sample  # doctest: +REMOTE_DATA

    >>> my_map = sunpy.map.Map(sunpy.data.sample.AIA_171_IMAGE)  # doctest: +REMOTE_DATA

The variable my_map is an `~sunpy.map.sources.AIAMap` object. To create one from a
local FITS file try the following::

    >>> my_map = sunpy.map.Map('/mydirectory/mymap.fits')   # doctest: +SKIP

SunPy should automatically detects the type of file (e.g. FITS), what instrument it is
associated with (e.g. AIA, EIT, LASCO) and will automatically look in the
appropriate places for the FITS keywords it needs to interpret the coordinate
system. If the type of FITS file is not recognized then SunPy will try some
default FITS keywords and return a `~sunpy.map.GenericMap` but results
may vary. SunPy can also create maps from the jpg2000 files from
`helioviewer.org <https://helioviewer.org/>`_.

Creating Custom Maps
====================
It is also possible to create maps using custom data (e.g. from a simulation or an observation
from a data source that is not explicitly supported in SunPy.) To do this you need to provide
`sunpy.map.Map` with both the data array as well as appropriate
meta information. The meta information is important as it informs the `sunpy.map.Map`
of the correct coordinate information associated with the data array. The meta information should be provided to
`sunpy.map.Map` in the form of a header as a `dict` or `~sunpy.util.MetaDict`.

The keys that are required for the header information follows the `FITS standard <https://fits.gsfc.nasa.gov/fits_dictionary.html>`_. SunPy now provides a map header helper function to assist the user in creating a header that contains the correct meta information
to generate a `sunpy.map.Map`.

The helper functionality includes a `~sunpy.map.meta_keywords` function
that will return a `dict` of all the current meta keywords and their descriptions currently used by
`sunpy.map.Map` to make a map::

    >>> from sunpy.map import meta_keywords

    >>> meta_keywords() # doctest: +SKIP
    {'cunit1': 'Units of the coordinate increments along naxis1 e.g. arcsec **required',
     'cunit2': 'Units of the coordinate increments along naxis2 e.g. arcsec **required',
     'crval1': 'Coordinate value at reference point on naxis1 **required'
     ...

There is also functionality also includes a utility function
`~sunpy.map.make_fitswcs_header` that will return a header with the
appropiate FITS keywords once the map data array and an `astropy.coordinates.SkyCoord` or `sunpy.coordinates.frames`
is passed. The `astropy.coordinates.SkyCoord` is defined by the user, and contains information on the reference frame,
reference coordinate and observer location. The function returns a `sunpy.util.MetaDict`.
The `astropy.coordinates.SkyCoord` or `sunpy.coordinates.frames` must contain an observation time.

The `~sunpy.map.make_fitswcs_header` function also takes optional keywords arguments including ``reference_pixel`` and ``scale`` which describe the pixel coordinate at the reference coordinate (defined by the `~astropy.coordinates.SkyCoord`) and the spatial scale of the pixels, respectively. If neither of these are given their values default to the center of the data array and 1 arcsec, respectively.

Here's an example of creating a header from some generic data and an `astropy.coordinates.SkyCoord`::


    >>> import numpy as np
    >>> import astropy.units as u
    >>> from sunpy.coordinates import frames
    >>> from astropy.coordinates import SkyCoord

    >>> data = np.arange(0,100).reshape(10,10)
    >>> coord = SkyCoord(0*u.arcsec, 0*u.arcsec, obstime = '2013-10-28', observer = 'earth', frame = frames.Helioprojective)
    >>> header = sunpy.map.make_fitswcs_header(data, coord)
    >>> for key, value in header.items():
    ...     print(f"{key}: {value}")
    wcsaxes: 2
    crpix1: 5.5
    crpix2: 5.5
    cdelt1: 1.0
    cdelt2: 1.0
    cunit1: arcsec
    cunit2: arcsec
    ctype1: HPLN-TAN
    ctype2: HPLT-TAN
    crval1: 0.0
    crval2: 0.0
    ...
    date-obs: 2013-10-28T00:00:00.000
    rsun_ref: 695700000.0
    dsun_obs: 148644585949.49
    hgln_obs: 0.0
    hglt_obs: 4.7711570596394


From this we can see now that the function returned a `sunpy.util.MetaDict` that populated
the standard FITS keywords with information provided by the passed `astropy.coordinates.SkyCoord`,
and the data array. Since the ``reference_pixel`` and keywords were not passed in the example above, the
values of ``crpix`` and ``cdelt`` were set to the default values.

These keywords can be passed to the function in the form of an `astropy.units.Quantity` with associated units.
Here's another example of passing ``reference_pixel`` and ``scale`` to the function::

    >>> header = sunpy.map.make_fitswcs_header(data, coord,
    ...                                        reference_pixel=u.Quantity([5, 5]*u.pixel),
    ...                                        scale=u.Quantity([2, 2] *u.arcsec/u.pixel))
    >>> for key, value in header.items():
    ...     print(f"{key}: {value}")
    wcsaxes: 2
    crpix1: 6.0
    crpix2: 6.0
    cdelt1: 2.0
    cdelt2: 2.0
    cunit1: arcsec
    cunit2: arcsec
    ctype1: HPLN-TAN
    ctype2: HPLT-TAN
    crval1: 0.0
    crval2: 0.0
    ...
    date-obs: 2013-10-28T00:00:00.000
    rsun_ref: 695700000.0
    dsun_obs: 148644585949.49
    hgln_obs: 0.0
    hglt_obs: 4.7711570596394

As we can see, a list of WCS and observer meta information is contained within the generated headers,
however we may want to include other meta information including the observatory name, the wavelength and
waveunit of the observation. Any of the keywords listed in ``header_helper.meta_keywords`` can be passed
to the `~sunpy.map.make_fitswcs_header` and will then populate the returned MetaDict header.
Furthermore, the following observation keywords can be passed to the `~sunpy.map.make_fitswcs_header`
function and will be translated to the FITS standard: ``observtory``, ``instrument``,``telescope``, ``wavelength``, ``exposure``.

An example of creating a header with these additional keywords::

    >>> header = sunpy.map.make_fitswcs_header(data, coord,
    ...                                        reference_pixel = u.Quantity([5, 5]*u.pixel),
    ...                                        scale = u.Quantity([2, 2] *u.arcsec/u.pixel),
    ...                                        telescope = 'Test case', instrument = 'UV detector',
    ...                                        wavelength = 1000*u.angstrom)
    >>> header  # doctest: +SKIP
    MetaDict([('wcsaxes', 2),
          ('crpix1', 5.0),
          ('crpix2', 5.0),
          ('cdelt1', <Quantity 2. arcsec2 / pix2>),
          ('cdelt2', <Quantity 2. arcsec2 / pix2>),
          ('cunit1', Unit("arcsec")),
          ('cunit2', Unit("arcsec")),
          ('ctype1', 'HPLN-TAN'),
          ('ctype2', 'HPLT-TAN'),
          ('crval1', 0.0),
          ('crval2', 0.0),
          ...
          ('date-obs', '2013-10-28T00:00:00.000'),
          ('hgln_obs', 0.0),
          ('hglt_obs', 4.7711570596394015),
          ('dsun_obs', 148644585949.4918),
          ('rsun_ref', 695700.0),
          ('rsun_obs', 965.3829548285768),
          ('instrume', 'Test case'),
          ('wavelnth', 1000),
          ('detector', 'UV detector'),
          ('waveunit', 'angstrom')])

From these header MetaDict's that are generated, we can now create a custom map::

    >>> my_map = sunpy.map.Map(data, header) # doctest: +SKIP
    >>> my_map.peek() # doctest: +SKIP

Inspecting maps
===============
A map contains a number of data-associated attributes. To get a quick look at
your map simply type::

    >>> my_map = sunpy.map.Map(sunpy.data.sample.AIA_171_IMAGE)  # doctest: +REMOTE_DATA
    >>> my_map  # doctest: +REMOTE_DATA
    <sunpy.map.sources.sdo.AIAMap object at ...>
    SunPy Map
    ---------
    Observatory:                 SDO
    Instrument:          AIA 3
    Detector:            AIA
    Measurement:                 171.0 Angstrom
    Wavelength:          171.0 Angstrom
    Observation Date:    2011-06-07 06:33:02
    Exposure Time:               0.234256 s
    Dimension:           [1024. 1024.] pix
    Coordinate System:   helioprojective
    Scale:                       [2.402792 2.402792] arcsec / pix
    Reference Pixel:     [511.5 511.5] pix
    Reference Coord:     [3.22309951 1.38578135] arcsec
    array([[ -95.92475  ,    7.076416 ,   -1.9656711, ..., -127.96519  ,
            -127.96519  , -127.96519  ],
           [ -96.97533  ,   -5.1167884,    0.       , ...,  -98.924576 ,
            -104.04137  , -127.919716 ],
           [ -93.99607  ,    1.0189276,   -4.0757103, ...,   -5.094638 ,
             -37.95505  , -127.87541  ],
           ...,
           [-128.01454  , -128.01454  , -128.01454  , ..., -128.01454  ,
            -128.01454  , -128.01454  ],
           [-127.899666 , -127.899666 , -127.899666 , ..., -127.899666 ,
            -127.899666 , -127.899666 ],
           [-128.03072  , -128.03072  , -128.03072  , ..., -128.03072  ,
            -128.03072  , -128.03072  ]], dtype=float32)

This will show a representation of the data as well as some of its associated
attributes. A number of other attributes are also available, for example the
`~sunpy.map.GenericMap.date`, `~sunpy.map.GenericMap.exposure_time`,
`~sunpy.map.GenericMap.center` and others (see `~sunpy.map.GenericMap`)::

    >>> map_date = my_map.date  # doctest: +REMOTE_DATA
    >>> map_exptime = my_map.exposure_time  # doctest: +REMOTE_DATA
    >>> map_center = my_map.center  # doctest: +REMOTE_DATA

To get a list of all of the attributes check the documentation by typing::

    >>> help(my_map)  # doctest: +SKIP

Many attributes and functions of the map classes accept and return
`~astropy.units.quantity.Quantity` or `~astropy.coordinates.SkyCoord` objects,
please refer to :ref:`units-coordinates-sunpy` for more details.

The meta data for the map is accessed by ::

    >>> header = my_map.meta  # doctest: +REMOTE_DATA

This references the meta data dictionary with the header information as read
from the source file.

Getting at the data
===================
The data in a SunPy Map object is accessible through the
`~sunpy.map.GenericMap.data` attribute.  The data is implemented as a
NumPy `~numpy.ndarray`, so for example, to get
the 0th element in the array ::

    >>> my_map.data[0, 0]  # doctest: +REMOTE_DATA
    -95.92475
    >>> my_map.data[0][0]  # doctest: +REMOTE_DATA
    -95.92475

One important fact to remember is that the first
index is for the y direction while the second index is for the x direction.
For more information about indexing please refer to the
`Numpy documentation <https://docs.scipy.org/doc/numpy-dev/user/quickstart.html#indexing-slicing-and-iterating>`_.

Data attributes like `~numpy.ndarray.dtype` and
`~sunpy.map.GenericMap.dimensions` are accessible through
the SunPyGenericMap object ::

    >>> my_map.dimensions  # doctest: +REMOTE_DATA
    PixelPair(x=<Quantity 1024. pix>, y=<Quantity 1024. pix>)
    >>> my_map.dtype  # doctest: +REMOTE_DATA
    dtype('float32')

Here the dimensions attribute is similar to the `~numpy.ndarray.shape`
attribute, however returning an `~astropy.units.quantity.Quantity`.

If you'd like to use the data in a SunPy `~sunpy.map.GenericMap` object
elsewhere, you can use either of the following::

    >>> var = my_map.data  # doctest: +REMOTE_DATA
    >>> var = my_map.data.copy()  # doctest: +REMOTE_DATA

Python makes use of pointers so if you want to alter the data and keep the
original data in the map intact make sure to copy it.

To create a complete copy of a Map object that is entirely independent of the original,
use the built-in `copy.deepcopy`
method, like so::

    >>> import copy   # doctest: +REMOTE_DATA
    >>> my_map_deepcopy = copy.deepcopy(my_map)   # doctest: +REMOTE_DATA

A deepcopy ensures that any changes in the original Map object are not reflected in the
copied object and vice versa. Note that this is different from simply copying the data of
the Map object - this copies all of the other attributes and methods as well.

Some basic statistical functions on the data array are also passed through to Map
objects::

    >>> my_map.min()  # doctest: +REMOTE_DATA
    -129.78036
    >>> my_map.max()  # doctest: +REMOTE_DATA
    192130.17
    >>> my_map.mean()  # doctest: +REMOTE_DATA
    427.02252

but you can also access all the other `~numpy.ndarray` functions and attributes
by accessing the data array directly. For example::

    >>> my_map.data.std()  # doctest: +REMOTE_DATA
    826.41016

Plotting
========
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

   If the `astropy.visualization.wcsaxes` package is not used (it is used by
   default) the `~sunpy.map.GenericMap.plot()` and
   `~sunpy.map.GenericMap.peek()` methods assume that the data is not rotated,
   i.e. the solar y axis is oriented with the columns of the array. If this
   condition is not met (in the metadata), when the map is plotted a warning
   will be issued. You can create an oriented map by using
   `~sunpy.map.GenericMap.rotate()` before you plot the Map.

Plotting Keywords
-----------------

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

You can view or make changes to the default settings through the ``sunpy.map.GenericMap.plot_settings``
dictionary. In the following example we change the title of the plot by changing the
``sunpy.map.GenericMap.plot_settings`` property.

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
---------------------------

Image data is generally shown in false color in order to better identify it or
to better visualize structures in the image. Matplotlib handles this colormapping
process through the `~matplotlib.colors` module. This process involves two steps:
the data array is first mapped onto the range 0-1 using an instance of
`~matplotlib.colors.Normalize` or a subclass; then this number is mapped to a
color using an instance of a subclass of a `~matplotlib.colors.Colormap`.

SunPy provides the colormaps for each mission as defined by the mission teams.
The Map object chooses the appropriate colormap for you when it is created as
long as it recognizes the instrument. To see what colormaps are available::

    >>> import sunpy.visualization.colormaps as cm
    >>> cm.cmlist.keys()
    dict_keys(['goes-rsuvi94', 'goes-rsuvi131', 'goes-rsuvi171', 'goes-rsuvi195',
    'goes-rsuvi284', 'goes-rsuvi304', 'sdoaia94', 'sdoaia131', 'sdoaia171',
    'sdoaia193', 'sdoaia211', 'sdoaia304', 'sdoaia335', 'sdoaia1600', 'sdoaia1700',
    'sdoaia4500', 'sohoeit171', 'sohoeit195', 'sohoeit284', 'sohoeit304', 'soholasco2',
    'soholasco3', 'sswidlsoholasco2', 'sswidlsoholasco3', 'stereocor1',
    'stereocor2', 'stereohi1', 'stereohi2', 'yohkohsxtal',
    'yohkohsxtwh', 'hinodexrt', 'hinodesotintensity', 'trace171', 'trace195',
    'trace284', 'trace1216', 'trace1550', 'trace1600', 'trace1700', 'traceWL',
    'hmimag', 'irissji1330', 'irissji1400', 'irissji1600', 'irissji2796',
    'irissji2832', 'irissji5000', 'irissjiFUV', 'irissjiNUV', 'irissjiSJI_NUV', 'kcor',
    'rhessi', 'std_gamma_2'])

The SunPy colormaps are registered with matplotlib so you can grab them like
you would any other colormap::

    >>> import matplotlib.pyplot as plt
    >>> import sunpy.visualization.colormaps

You need to import `sunpy.visualization.colormaps` or `sunpy.map` for this to work::

    >>> cmap = plt.get_cmap('sdoaia171')


The following plot shows off all of the colormaps.

.. plot::

    import matplotlib.pyplot as plt
    import sunpy.visualization.colormaps as cm
    cm.show_colormaps()

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
    smap.plot(cmap=cmap)
    plt.colorbar()
    plt.show()

or you can just change the colormap for the map itself as follows::

    >>> smap.plot_settings['cmap'] = plt.get_cmap('sohoeit171')  # doctest: +SKIP

The normalization is also set automatically and is chosen so that all the
data from minimum to maximum is displayed as best as possible for most cases.
This means that it is never necessary to touch the data such as applying a function
such sqrt or log to the data to make your plot look good.
There are many normalizations available from matplotlib such as `~matplotlib.colors.LogNorm`. Other
`more exotic normalizations <https://docs.astropy.org/en/stable/visualization/index.html>`_ are also
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

    fig = plt.figure(figsize=(4, 9))

    ax1 = fig.add_subplot(2, 1, 1, projection=smap)
    smap.plot(norm=colors.Normalize(), title='Linear normalization')
    plt.colorbar()

    ax2 = fig.add_subplot(2, 1, 2, projection=smap)
    smap.plot(norm=colors.LogNorm(), title='Logarithmic normalization')
    plt.colorbar()

    plt.show()

Note how the color in the colorbar does not change since these two maps share
the same colormap while the data values associated with each color do because
the normalization is different.

Masking and Clipping Data
=========================
It is often necessary for the purposes of display or otherwise to ignore certain
data in an image. For example, a large data value could be due to
cosmic ray hits and should be ignored. The most straightforward way to ignore
this kind of data in plots without altering the data is to clip it. This can be achieved
very easily by using the ``clip_interval`` keyword. For example::

    >>> import astropy.units as u
    >>> smap.plot(clip_interval=(1, 99.5)*u.percent)  #doctest: +SKIP

This clips out the dimmest 1% of pixels and the brightest 0.5% of pixels.  With those outlier
pixels clipped, the resulting image makes better use of the full range of colors.
If you'd like to see what areas of your images got clipped, you can modify the colormap::

    >>> cmap = map.cmap  # doctest: +SKIP
    >>> cmap.set_over('blue')  # doctest: +SKIP
    >>> cmap.set_under('green')  # doctest: +SKIP

This will color the areas above and below in red and green respectively
(similar to this `example <https://matplotlib.org/examples/pylab_examples/image_masked.html>`_).
You can use the following colorbar command to display these choices::

    >>> plt.colorbar(extend='both')   # doctest: +SKIP

Here is an example of this put to use on an AIA image.

.. plot::
    :include-source:

    import astropy.units as u
    import matplotlib.pyplot as plt

    import sunpy.map
    import sunpy.data.sample

    smap = sunpy.map.Map(sunpy.data.sample.AIA_171_IMAGE)
    cmap = smap.cmap
    cmap.set_over('blue')
    cmap.set_under('green')

    fig = plt.figure(figsize=(12, 4))

    ax1 = fig.add_subplot(1, 2, 1, projection=smap)
    smap.plot(title='Without clipping')
    plt.colorbar()

    ax2 = fig.add_subplot(1, 2, 2, projection=smap)
    smap.plot(clip_interval=(1, 99.5)*u.percent, title='With clipping')
    plt.colorbar(extend='both')

    plt.show()

Another approach to clipping data is to specify explicit values for the minimum and maximum pixel
values using the plotting keywords ``vmin`` and ``vmax``.

Clipping excludes data that has extreme values, but there can be other forms of bad data.
A mask is a boolean
array and so can give you much more fine-grained control over what is not being
displayed.  A `~numpy.ma.MaskedArray`
is a subclass of a numpy array so it has all of the same properties with the
addition of an associated boolean array which holds the mask.
See `this example <https://docs.sunpy.org/en/stable/generated/gallery/computer_vision_techniques/mask_disk.html>`_ in our gallery.

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
    txt = r"min={min}, max={max}, $\mu$={mean}, $\sigma$={std}".format(min=int(smap.min()),
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
    txt = r"min={min}, max={max}, $\mu$={mean}, $\sigma$={std}".format(min=int(smap.min()),
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
    txt = r"min={min}, max={max}, $\mu$={mean}, $\sigma$={std}".format(min=int(smap.min()),
                                                                       max=int(smap.max()),
                                                                       mean=int(smap.mean()),
                                                                       std=int(smap.std()))
    plt.text(-600, 1500, txt, color='white')
    norm = colors.Normalize()
    smap.plot(norm = norm)
    plt.colorbar(extend='both')


Composite Maps and Overlaying Maps
==================================

The `~sunpy.map.Map` method described above can also handle a list of maps. If a series of maps
are supplied as inputs, `~sunpy.map.Map` will return a list of maps as the output.  However,
if the 'composite' keyword is set to True, then a `~sunpy.map.CompositeMap` object is
returned.  This is useful if the maps are of a different type (e.g. different
instruments).  For example, to create a simple composite map::

    >>> my_maps = sunpy.map.Map(sunpy.data.sample.EIT_195_IMAGE, sunpy.data.sample.RHESSI_IMAGE, composite=True)  # doctest: +REMOTE_DATA

A `~sunpy.map.CompositeMap` is different from a regular SunPy `~sunpy.map.GenericMap` object and therefore
different associated methods. To list which maps are part of your composite map use::

    >>> my_maps.list_maps()  # doctest: +REMOTE_DATA
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
=====================
Part of the philosophy of the map object is to provide most of the basic
functionality that a scientist would want therefore a map also contains a number
of map-specific methods such as resizing a map or grabbing a subview. To get
a list of the methods available for a map type::

    >>> help(my_map)  # doctest: +SKIP

and check out the methods section!

MapSequences
============
A `~sunpy.map.MapSequence` is an ordered list of maps.  By default, the maps are ordered by
their observation date, from earlier maps to later maps. A `~sunpy.map.MapSequence` can be
created by supplying multiple existing maps::

    >>> map1 = sunpy.map.Map(sunpy.data.sample.AIA_171_IMAGE)  # doctest: +REMOTE_DATA
    >>> map2 = sunpy.map.Map(sunpy.data.sample.EIT_195_IMAGE)  # doctest: +REMOTE_DATA
    >>> mc = sunpy.map.Map([map1, map2], sequence=True)  # doctest: +REMOTE_DATA

or by providing a directory full of image files::

    >>> mc = sunpy.map.Map('path/to/my/files/*.fits', sequence=True)   #  doctest: +SKIP

The earliest map in the MapSequence can be accessed by simply indexing the maps
list::

    >>> mc.maps[0]   # doctest: +SKIP

MapSequences can hold maps that have different shapes.  To test if all the
maps in a `~sunpy.map.MapSequence` have the same shape::

    >>> mc.all_maps_same_shape()  # doctest: +REMOTE_DATA
    True

It is often useful to return the image data in a `~sunpy.map.MapSequence` as a single
three dimensional Numpy `~numpy.ndarray`::

    >>> mc.as_array()   # doctest: +SKIP

Note that an array is returned only if all the maps have the same
shape.  If this is not true, an error (ValueError) is returned.  If all the
maps have nx pixels in the x-direction, and ny pixels in the y-direction,
and there are n maps in the MapSequence, the `~numpy.ndarray` array that is
returned has shape (ny, nx, n).  The data of the first map in the `~sunpy.map.MapSequence`
appears in the `~numpy.ndarray` in position ``[:, :, 0]``, the data of second map in
position ``[:, :, 1]``, and so on.  The order of maps in the `~sunpy.map.MapSequence` is
reproduced in the returned `~numpy.ndarray`.

The meta data from each map can be obtained using::

    >>> mc.all_meta()   # doctest: +SKIP

This returns a list of map meta objects that have the same order as
the maps in the `~sunpy.map.MapSequence`.

Coalignment of MapSequences
===========================
A typical data preparation step when dealing with time series of images is to
coalign images taken at different times so that features in different images
remain in the same place.  A common approach to this problem is
to take a representative template that contains the features you are interested
in, and match that to your images.  The location of the best match tells you
where the template is in your image.  The images are then shifted to the
location of the best match.  This aligns your images to the position of the
features in your representative template.

SunPy provides a function to coalign the maps inside the `~sunpy.map.MapSequence`.
The implementation of this functionality requires the installation of the
scikit-image library, a commonly used image processing library.
To coalign a `~sunpy.map.MapSequence`, simply import
the function and apply it to your `~sunpy.map.MapSequence`::

    >>> from sunpy.image.coalignment import mapsequence_coalign_by_match_template
    >>> coaligned = mapsequence_coalign_by_match_template(mc)  # doctest: +REMOTE_DATA

This will return a new `~sunpy.map.MapSequence`, coaligned to a template extracted from the
center of the first map in the `~sunpy.map.MapSequence`, with the map dimensions clipped as
required.  The coalignment algorithm provides many more options for handling
the coalignment of `~sunpy.map.MapSequence` type::

    >>> help(mapsequence_coalign_by_match_template)   # doctest: +SKIP

for a full list of options and functionality.

If you just want to calculate the shifts required to compensate for solar
rotation relative to the first map in the `~sunpy.map.MapSequence` without applying them, use::

    >>> from sunpy.image.coalignment import calculate_match_template_shift
    >>> shifts = calculate_match_template_shift(mc)  # doctest: +REMOTE_DATA

This is the function used to calculate the shifts in `~sunpy.map.MapSequence` coalignment
function above.  Please see `~sunpy.image.coalignment.calculate_match_template_shift` to learn more about its features.
Shifts calculated using calculate_match_template_shift can be passed directly
to the coalignment function.


Compensating for solar rotation in MapSequences
===============================================
Often a set of solar image data consists of fixing the pointing of a
field of view for some time and observing.  Features on the Sun will
rotate according to the Sun's rotation.

A typical data preparation step when dealing with time series of these
types of images is to shift the images so that features do not appear
to move across the field of view.  This requires taking in to account
the rotation of the Sun.  The Sun rotates differentially, depending on
latitude, with features at the equator moving faster than features at
the poles.

SunPy provides a function to shift images in `~sunpy.map.MapSequence` following solar
rotation.  This function shifts an image according to the solar
differential rotation calculated at the latitude of the center of the
field of view.  The image is not *differentially* rotated.  This
function is useful for de-rotating images when the effects of
differential rotation in the `~sunpy.map.MapSequence` can be ignored (for example, if
the spatial extent of the image is small, or when the duration of the
`~sunpy.map.MapSequence` is small; deciding on what 'small' means depends on your
application).

To apply this form of solar derotation to a `~sunpy.map.MapSequence`, simply import the
function and apply it to your `~sunpy.map.MapSequence`::

    >>> from sunpy.physics.solar_rotation import mapsequence_solar_derotate
    >>> derotated = mapsequence_solar_derotate(mc)  # doctest: +SKIP

For more info see `~sunpy.physics.solar_rotation.mapsequence_solar_derotate`.

If you just want to calculate the shifts required to compensate for solar
rotation relative to the first map in the `~sunpy.map.MapSequence` without applying them, use::

    >>> from sunpy.physics.solar_rotation import calculate_solar_rotate_shift
    >>> shifts = calculate_solar_rotate_shift(mc)  # doctest: +SKIP

Please consult the docstring of the `~sunpy.image.coalignment.mapsequence_coalign_by_match_template` function in order to learn about the features of this function.
