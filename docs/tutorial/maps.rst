.. _map_guide:

****
Maps
****

In this section of the tutorial, you will learn about the `Map <sunpy.map.GenericMap>` object, the primary data structure provided by sunpy.
Map objects hold two-dimensional data along with the accompanying metadata.
They can be used with any two-dimensional data array with two spatial axes and FITS-compliant metadata.
As you will see in this tutorial, the primary advantage of using a Map object over just a bare NumPy array is the ability to perform coordinate aware operations on the underlying array, such as rotating the Map to remove the roll angle of an instrument or cropping a Map to a specific field of view.
Additionally, because Map is capable of ingesting data from many different sources, it provides a common interface for any two-dimensional data product.

By the end of this tutorial, you will learn how to create a Map, access the underlying data and metadata, and the basics visualizing a Map.
Additionally, you will learn how to combine multiple Maps into a `~sunpy.map.MapSequence` or a `~sunpy.map.CompositeMap`.
Lastly, you will learn how to create Maps from data sources not supported in sunpy.

.. note::

    In this section and in :ref:`timeseries_guide`, we will use the sample data included with sunpy.
    These data are primarily useful for demonstration purposes or simple debugging.
    These files have names like ``sunpy.data.sample.AIA_171_IMAGE`` and ``sunpy.data.sample.RHESSI_IMAGE`` and are automatically downloaded to your computer as you need them.
    Once downloaded, these sample data files will be paths to their location on your computer.

.. _creating-maps:

Creating Maps
=============

To create a `~sunpy.map.Map` from a sample AIA image,

.. code-block:: python

    >>> import sunpy.map
    >>> import sunpy.data.sample  # doctest: +REMOTE_DATA
    >>> my_map = sunpy.map.Map(sunpy.data.sample.AIA_171_IMAGE)  # doctest: +REMOTE_DATA

sunpy automatically detects the type of file (e.g. FITS) as well as the the instrument associated with it (e.g. AIA, EIT, LASCO).
To make sure this has all worked correctly, we can take a quick look at ``my_map``,

.. code-block:: python

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

This should show a representation of the data as well as some of its associated attributes.
If you are in a Jupyter Notebook, this will show a rich HTML view of the table along with several quick-look plots.
Otherwise, you can use the :meth:`~sunpy.map.GenericMap.quicklook` method to see these quick-look plots.

.. code-block:: python

    >>> my_map.quicklook()  # doctest: +SKIP

.. generate:: html
    :html_border:

    import sunpy.map
    import sunpy.data.sample
    my_map = sunpy.map.Map(sunpy.data.sample.AIA_171_IMAGE)
    print(my_map._repr_html_())

.. _inspecting-maps:

Inspecting Map Metadata
=======================

The metadata for a Map is exposed via attributes on the Map.
These attributes can be accessed by typing `my_map.<attribute-name>`.
For example, to access the date of the observation,

.. code-block:: python

    >>> my_map.date  # doctest: +REMOTE_DATA
    <Time object: scale='utc' format='isot' value=2011-06-07T06:33:02.770>

Notice that this is an `~astropy.time.Time` object which we discussed in the previous section :ref:`time-in-sunpy`.
Similarly, we can access the exposure time of the image,

.. code-block:: python

    >>> my_map.exposure_time  # doctest: +REMOTE_DATA
    <Quantity 0.234256 s>

Notice that this returns an `~astropy.units.Quantity` object which we discussed in the previous section :ref:`units-sunpy`.
The full list of attributes can be found on `~sunpy.map.GenericMap`.
These metadata attributes are all derived from the underlying FITS metadata, but are represented as rich Python objects, rather than simple strings or floating point numbers.

.. _map-data:

Map Data
========

The data in a Map is stored as a `numpy.ndarray` object and is accessible through the `~sunpy.map.GenericMap.data` attribute:

.. code-block:: python

    >>> my_map.data  # doctest: +REMOTE_DATA
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

This array can then be indexed like any other NumPy array.
For example, to get the 0th element in the array:

.. code-block:: python

    >>> my_map.data[0, 0]  # doctest: +REMOTE_DATA
    -95.92475

The first index corresponds to the y direction and the second to the x direction in the two-dimensional pixel coordinate system.
For more information about indexing, please refer to the `numpy documentation <https://numpy.org/doc/stable/user/basics.indexing.html#indexing-on-ndarrays>`__.

Data attributes like dimensionality and type are also accessible as attributes on ``my_map``:

.. code-block:: python

    >>> my_map.dimensions  # doctest: +REMOTE_DATA
    PixelPair(x=<Quantity 1024. pix>, y=<Quantity 1024. pix>)
    >>> my_map.dtype  # doctest: +REMOTE_DATA
    dtype('float32')

Additional, there are several methods that provide basic summary statistics of the data:

.. code-block:: python

    >>> my_map.min()  # doctest: +REMOTE_DATA
    -129.78036
    >>> my_map.max()  # doctest: +REMOTE_DATA
    192130.17
    >>> my_map.mean()  # doctest: +REMOTE_DATA
    427.02252

.. _coordinates-wcs-maps:

Maps, Coordinates, and the World Coordinate System
==================================================

In :ref:`coordinates-sunpy`, you learned how to define coordinates with `~astropy.coordinates.SkyCoord` using different solar coordinate frames.
The coordinate frame of a Map is provided as an attribute,

.. code-block:: python

    >>> my_map.coordinate_frame  # doctest: +REMOTE_DATA
    <Helioprojective Frame (obstime=2011-06-07T06:33:02.770, rsun=696000.0 km, observer=<HeliographicStonyhurst Coordinate (obstime=2011-06-07T06:33:02.770, rsun=696000.0 km): (lon, lat, radius) in (deg, deg, m)
        (-0.00406308, 0.04787238, 1.51846026e+11)>)>

This tells us that the coordinate system of the image is Helioprojective (HPC) and that it is defined by an observer at a particular location.
This observer coordinate is also provided as an attribute,

.. code-block:: python

    >>> my_map.observer_coordinate  # doctest: +REMOTE_DATA
    <SkyCoord (HeliographicStonyhurst: obstime=2011-06-07T06:33:02.770, rsun=696000.0 km): (lon, lat, radius) in (deg, deg, m)
        (-0.00406308, 0.04787238, 1.51846026e+11)>

This tells us the location of the spacecraft, in this case SDO, when it recorded this partiuclar observation, as derived from the FITS metadata.
Map has several additional coordinate-related attributes that provide the coordinates of the center and corners of the Map,

.. code-block:: python

    >>> my_map.center  # doctest: +REMOTE_DATA
    <SkyCoord (Helioprojective: obstime=2011-06-07T06:33:02.770, rsun=696000.0 km, observer=<HeliographicStonyhurst Coordinate (obstime=2011-06-07T06:33:02.770, rsun=696000.0 km): (lon, lat, radius) in (deg, deg, m)
        (-0.00406308, 0.04787238, 1.51846026e+11)>): (Tx, Ty) in arcsec
        (3.22309951, 1.38578135)>
    >>> my_map.bottom_left_coord  # doctest: +REMOTE_DATA
    <SkyCoord (Helioprojective: obstime=2011-06-07T06:33:02.770, rsun=696000.0 km, observer=<HeliographicStonyhurst Coordinate (obstime=2011-06-07T06:33:02.770, rsun=696000.0 km): (lon, lat, radius) in (deg, deg, m)
        (-0.00406308, 0.04787238, 1.51846026e+11)>): (Tx, Ty) in arcsec
        (-1228.76466158, -1224.62447509)>
    >>> my_map.top_right_coord  # doctest: +REMOTE_DATA
    <SkyCoord (Helioprojective: obstime=2011-06-07T06:33:02.770, rsun=696000.0 km, observer=<HeliographicStonyhurst Coordinate (obstime=2011-06-07T06:33:02.770, rsun=696000.0 km): (lon, lat, radius) in (deg, deg, m)
        (-0.00406308, 0.04787238, 1.51846026e+11)>): (Tx, Ty) in arcsec
        (1235.21095899, 1227.39598836)>

But what if we wanted to know what pixel these physical coordinates correspond to?
Thankfully, each Map has an associated World Coordinate System, or WCS, which is derived from the underlying metadata and expressed as an `astropy.wcs.WCS` object.
The WCS is accessible as an attribute:

.. code-block:: python

    >>> my_map.wcs  # doctest: +REMOTE_DATA
    WCS Keywords
    <BLANKLINE>
    Number of WCS axes: 2
    CTYPE : 'HPLN-TAN'  'HPLT-TAN'
    CRVAL : 0.00089530541880571  0.00038493926472939
    CRPIX : 512.5  512.5
    PC1_1 PC1_2  : 0.99999706448085  0.0024230207763071
    PC2_1 PC2_2  : -0.0024230207763071  0.99999706448085
    CDELT : 0.00066744222222222  0.00066744222222222
    NAXIS : 1024  1024

WCS is a fairly complex topic, but all we need to know for now is that the WCS provides the transformation between the pixel coordinates of the image and physical or "world" coordinates.
In particular, we will only focus on two methods: `~astropy.wcs.WCS.world_to_pixel` and `~astropy.wcs.WCS.pixel_to_world`.
First, let's find the pixel location corresponding to the center of the Map,

.. code-block:: python

    >>> center_pixel = my_map.wcs.world_to_pixel(my_map.center)  # doctest: +REMOTE_DATA
    >>> center_pixel  # doctest: +REMOTE_DATA
    (array(511.5), array(511.5))

Notice that these coordinates are not necessarily integers.
The corresponding pixel-to-world transformation should then give us back our center coordinate from above,

.. code-block:: python

    >>> my_map.wcs.pixel_to_world(center_pixel[0], center_pixel[1])  # doctest: +REMOTE_DATA
    <SkyCoord (Helioprojective: obstime=2011-06-07T06:33:02.770, rsun=696000.0 km, observer=<HeliographicStonyhurst Coordinate (obstime=2011-06-07T06:33:02.770, rsun=696000.0 km): (lon, lat, radius) in (deg, deg, m)
        (-0.00406308, 0.04787238, 1.51846026e+11)>): (Tx, Ty) in arcsec
        (3.22309951, 1.38578135)>

As another example, if we transform the center of the lower-left pixel to a world coordinate, it should correspond to bottom left coordinate from above,

.. code-block:: python

    >>> my_map.wcs.pixel_to_world(0, 0)  # doctest: +REMOTE_DATA
    <SkyCoord (Helioprojective: obstime=2011-06-07T06:33:02.770, rsun=696000.0 km, observer=<HeliographicStonyhurst Coordinate (obstime=2011-06-07T06:33:02.770, rsun=696000.0 km): (lon, lat, radius) in (deg, deg, m)
        (-0.00406308, 0.04787238, 1.51846026e+11)>): (Tx, Ty) in arcsec
        (-1228.76466158, -1224.62447509)>

These two methods are extremely useful when trying to understand which pixels correspond to which physical coordinates or when trying to locate the same physical location in images taken by separate spacecraft.

.. _plotting-maps:

Visualizing Maps
================

.. plot::
    :nofigs:
    :context: close-figs

    # This is here to put my_map in the scope of the plot directives.
    # This avoids repeating code in the example source code that is actually displayed.
    import sunpy.map
    import sunpy.data.sample
    my_map = sunpy.map.Map(sunpy.data.sample.AIA_171_IMAGE)

The Map object has a built-in plot method such that it is easy to quickly view your map.
sunpy makes use of `Matplotlib <https://matplotlib.org/>`_ for all of its plotting - as such, it tries to follow the Matplotlib plotting philosophy.
We refer the reader to the `Matplotlib usage documentation <https://matplotlib.org/stable/users/explain/api_interfaces.html>`__ to learn more about the Matplotlib and to become familiar with the basics.
To be consistent with Matplotlib, sunpy has developed a standard plotting interface which supports both simple and advanced Matplotlib usage.

Basic Plotting
--------------

For more advanced plotting the base sunpy objects also provide a `~sunpy.map.mapbase.GenericMap.plot` command.
This command is similar to the pyplot `~matplotlib.pyplot.imshow` command in that it will create a figure and axes object for you if you haven't already.

When you create a plot with `~sunpy.map.GenericMap.peek` or
`~sunpy.map.GenericMap.plot`, sunpy will use `astropy.visualization.wcsaxes` to
represent coordinates on the image accurately, for more information see
:ref:`wcsaxes-plotting`.

Using `~sunpy.map.GenericMap.plot` it is possible to customise the look of the
plot by combining sunpy and matplotlib commands, for example you can over plot
contours on the Map:

.. plot::
    :include-source:

    import matplotlib.pyplot as plt
    import astropy.units as u

    import sunpy.map
    import sunpy.data.sample

    aia_map = sunpy.map.Map(sunpy.data.sample.AIA_171_IMAGE)
    aia_map.plot()
    aia_map.draw_limb()

    # let's add contours as well
    aia_map.draw_contours([10,20,30,40,50,60,70,80,90] * u.percent)

    plt.colorbar()
    plt.show()


Plotting Keywords
-----------------

For Map plotting, `~matplotlib.pyplot.imshow` does most of the heavy lifting in the background while **sunpy** makes a number of choices for you (e.g. colortable, plot title).
Changing these defaults is made possible through two simple interfaces.
You can pass any `~matplotlib.pyplot.imshow` keyword into the plot command to override the defaults for that particular plot.
For example, the following plot changes the default colormap to use an inverse Grey color table.

.. plot::
    :include-source:

    import sunpy.map
    import sunpy.data.sample
    import matplotlib.pyplot as plt
    aia_map = sunpy.map.Map(sunpy.data.sample.AIA_171_IMAGE)
    fig = plt.figure()
    aia_map.plot(cmap=plt.cm.Greys_r)
    plt.colorbar()
    plt.show()

You can also view or make changes to the default settings through the ``sunpy.map.GenericMap.plot_settings`` dictionary.
See :ref:`sphx_glr_generated_gallery_plotting_map_editcolormap.py` for an example of this workflow for changing plot settings.


Changing the Colormap and Normalization
----------------------------------------

Image data is generally shown in false color in order to better identify it or to better visualize structures in the image.
Matplotlib handles this colormapping process through the `~matplotlib.colors` module.
First, the data array is mapped onto the range 0-1 using an instance of `~matplotlib.colors.Normalize` or a subclass.
Then, the data is mapped to a color using a `~matplotlib.colors.Colormap`.

**sunpy** provides colormaps for each mission as defined by the mission teams.
The Map object chooses the appropriate colormap for you when it is created as long as it recognizes the instrument.
To see what colormaps are available:

.. code-block:: python

    >>> import sunpy.visualization.colormaps as cm
    >>> cm.cmlist.keys()
    dict_keys(['goes-rsuvi94', 'goes-rsuvi131', 'goes-rsuvi171', 'goes-rsuvi195',
    'goes-rsuvi284', 'goes-rsuvi304', 'sdoaia94', 'sdoaia131', 'sdoaia171',
    ...

The **sunpy** colormaps are registered with Matplotlib so you can grab them like you would any other colormap:

.. code-block:: python

    >>> import matplotlib.pyplot as plt
    >>> import sunpy.visualization.colormaps
    >>> cmap = plt.get_cmap('sdoaia171')

See `~sunpy.visualization.colormaps` for a plot of all available colormaps.

If you want to override the built-in colormap, consider the following example which plots an AIA map using an EIT colormap.

.. plot::
    :include-source:

    import sunpy.map
    import sunpy.data.sample
    import matplotlib.pyplot as plt

    aia_map = sunpy.map.Map(sunpy.data.sample.AIA_171_IMAGE)
    cmap = plt.get_cmap('sohoeit171')

    fig = plt.figure()
    aia_map.plot(cmap=cmap)
    plt.colorbar()
    plt.show()

You can also change the colormap for the Map itself:

.. code-block:: python

    >>> smap.plot_settings['cmap'] = plt.get_cmap('sohoeit171')  # doctest: +SKIP

The normalization is set automatically so that all the data from minimum to maximum is displayed as best as possible.
Just like the colormap, the default normalization can be changed through the ``plot_settings`` dictionary or directly for the individual plot by passing a keyword argument.

Alternate normalizations are available from `matplotlib <https://matplotlib.org/stable/tutorials/colors/colormapnorms.html>`__ and `astropy <https://docs.astropy.org/en/stable/visualization/normalization.html>`__.
The following example shows the difference between a linear and logarithmic normalization on an AIA image.

.. plot::
    :include-source:

    import sunpy.map
    import sunpy.data.sample
    import matplotlib.pyplot as plt
    import matplotlib.colors as colors

    aia_map = sunpy.map.Map(sunpy.data.sample.AIA_171_IMAGE)

    fig = plt.figure(figsize=(4, 9))

    ax1 = fig.add_subplot(2, 1, 1, projection=aia_map)
    aia_map.plot(norm=colors.Normalize(), title='Linear normalization')
    plt.colorbar()

    ax2 = fig.add_subplot(2, 1, 2, projection=aia_map)
    aia_map.plot(norm=colors.LogNorm(), title='Logarithmic normalization')
    plt.colorbar()

    plt.show()

Note how the colorbar does not change since these two plots share the same colormap.
Meanwhile, the data values associated with each color do change because the normalization is different.


.. _wcsaxes-plotting:

Overlaying Coordinates
----------------------

By default :ref:`map` uses the `astropy.visualization.wcsaxes` module to improve
the representation of world coordinates, and calling
`~sunpy.map.GenericMap.plot` or `~sunpy.map.GenericMap.peek()` will use wcsaxes
for plotting. Unless a standard `matplotlib.axes.Axes` object is explicitly
created.

To explicitly create a `~astropy.visualization.wcsaxes.WCSAxes` instance do the
following ::

    >>> fig = plt.figure()   # doctest: +SKIP
    >>> ax = plt.subplot(projection=smap)   # doctest: +SKIP

when plotting on an `~astropy.visualization.wcsaxes.WCSAxes` axes, it will by
default plot in pixel coordinates, you can override this behavior and plot in
'world' coordinates by getting the transformation from the axes with
``ax.get_transform('world')``.

.. note::

    World coordinates are always in **degrees** so you will have to convert to degrees.

.. code-block:: python

    >>> aia_map.plot()   # doctest: +SKIP
    >>> ax.plot((100*u.arcsec).to_value(u.deg), (500*u.arcsec).to_value(u.deg),
    ...         transform=ax.get_transform('world'))   # doctest: +SKIP


In this next example, the `~matplotlib.figure.Figure` and
`~astropy.visualization.wcsaxes.WCSAxes` instances are created explicitly, and
then used to modify the plot.

Here we can plot a sunpy map, and also overplot some points defined in arcseconds, highlighting the advantage of using WCSAxes.

.. plot::
    :include-source:

    import matplotlib.pyplot as plt
    import astropy.units as u
    from astropy.coordinates import SkyCoord

    import sunpy.map
    import sunpy.data.sample

    aia_map = sunpy.map.Map(sunpy.data.sample.AIA_171_IMAGE)

    fig = plt.figure()
    # Provide the Map as a projection, which creates a WCSAxes object
    ax = plt.subplot(projection=aia_map)

    im = aia_map.plot()

    # Prevent the image from being re-scaled while overplotting.
    ax.set_autoscale_on(False)

    xc = [0,100,1000] * u.arcsec
    yc = [0,100,1000] * u.arcsec

    coords = SkyCoord(xc, yc, frame=aia_map.coordinate_frame)

    p = ax.plot_coord(coords, 'o')

    # Set title.
    ax.set_title('Custom plot with WCSAxes')

    plt.colorbar()
    plt.show()

It is possible to create the same plot, explicitly not using `~astropy.visualization.wcsaxes`, however, this will not have the features of `~astropy.visualization.wcsaxes` which include correct representation of rotation and plotting in different coordinate systems.


Check out the following foundational examples in the Example Gallery for plotting Maps:

* :ref:`sphx_glr_generated_gallery_plotting_aia_example.py`

* :ref:`sphx_glr_generated_gallery_map_plot_frameless_image.py`

* :ref:`sphx_glr_generated_gallery_plotting_wcsaxes_plotting_example.py`

* :ref:`sphx_glr_generated_gallery_plotting_map_editcolormap.py`

* :ref:`sphx_glr_generated_gallery_plotting_grid_plotting.py`


Clipping and Masking Data
-------------------------

It is often necessary to ignore certain data in an image.
For example, a large data value could be due to cosmic ray hits and should be ignored.
The most straightforward way to ignore this kind of data in plots, without altering the data, is to clip it.
This can be achieved very easily by using the ``clip_interval`` keyword. For example:

.. code-block:: python

    >>> import astropy.units as u
    >>> smap.plot(clip_interval=(1, 99.5)*u.percent)  #doctest: +SKIP

This clips out the dimmest 1% of pixels and the brightest 0.5% of pixels.
With those outlier pixels clipped, the resulting image makes better use of the full range of colors.
If you'd like to see what areas of your images got clipped, you can modify the colormap:

.. code-block:: python

    >>> cmap = map.cmap  # doctest: +SKIP
    >>> cmap.set_over('blue')  # doctest: +SKIP
    >>> cmap.set_under('green')  # doctest: +SKIP

This will color the areas above and below in red and green respectively (similar to this `matplotlib example <https://matplotlib.org/examples/pylab_examples/image_masked.html>`__).
You can use the following colorbar command to display these choices:

.. code-block:: python

    >>> plt.colorbar(extend='both')   # doctest: +SKIP

Here is an example of this put to use on an AIA image.

.. plot::
    :include-source:

    import astropy.units as u
    import matplotlib.pyplot as plt

    import sunpy.map
    import sunpy.data.sample

    smap = sunpy.map.Map(sunpy.data.sample.AIA_171_IMAGE)
    cmap = smap.cmap.copy()
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


Masking is another approach to ignoring certain data.
A mask is a boolean array that can give you fine-grained control over what is not being displayed.
The `~numpy.ma.MaskedArray` is a subclass of a NumPy array with the addition of an associated boolean array which holds the mask.
See the following two examples for applications of this technique:

* :ref:`sphx_glr_generated_gallery_computer_vision_techniques_mask_disk.py`

* :ref:`sphx_glr_generated_gallery_computer_vision_techniques_finding_masking_bright_pixels.py`

.. _map-sequences:

Map Sequences
=============

A `~sunpy.map.MapSequence` is an ordered list of maps.
By default, the maps are ordered by their observation date, from earliest to latest date.
A `~sunpy.map.MapSequence` can be created by supplying multiple existing maps:

.. code-block:: python

    >>> map1 = sunpy.map.Map(sunpy.data.sample.AIA_171_IMAGE)  # doctest: +REMOTE_DATA
    >>> map2 = sunpy.map.Map(sunpy.data.sample.EIT_195_IMAGE)  # doctest: +REMOTE_DATA
    >>> mc = sunpy.map.Map([map1, map2], sequence=True)  # doctest: +REMOTE_DATA

The earliest map in the MapSequence can be accessed by indexing the maps list:

.. code-block:: python

    >>> mc.maps[0]   # doctest: +SKIP

MapSequences can hold maps that have different shapes.
To test if all the maps in a `~sunpy.map.MapSequence` have the same shape:

.. code-block:: python

    >>> mc.all_maps_same_shape()  # doctest: +REMOTE_DATA
    True

It is often useful to return the image data in a `~sunpy.map.MapSequence` as a single three dimensional NumPy `~numpy.ndarray`:

.. code-block:: python

    >>> mc_array = mc.as_array()   # doctest: +REMOTE_DATA

Note that an array is returned only if all the maps have the same shape.
If this is not true, a `ValueError` is raised.
If all the maps have nx pixels in the x-direction, and ny pixels in the y-direction, and there are n maps in the MapSequence, the returned `~numpy.ndarray` array has shape (ny, nx, n).
The data of the first map in the `~sunpy.map.MapSequence` appears in the `~numpy.ndarray` in position ``[:, :, 0]``, the data of second map in position ``[:, :, 1]``, and so on.
The order of maps in the `~sunpy.map.MapSequence` is reproduced in the returned `~numpy.ndarray`.

The metadata from each map can be obtained using:

.. code-block:: python

    >>> mc.all_meta()   # doctest: +SKIP

This returns a list of map meta objects that have the same order as the maps in the `~sunpy.map.MapSequence`.

For information on coaligning images and compensating for solar rotation in Map Sequences, see the `sunkit-image example gallery <https://docs.sunpy.org/projects/sunkit-image/en/stable/generated/gallery/index.html>`__ and the `sunkit_image.coalignment` module.

.. _composite-maps:

Composite Maps and Overlaying Maps
==================================

The `~sunpy.map.Map` method can also handle a list of maps.
If a series of maps are supplied as inputs, `~sunpy.map.Map` will return a list of maps as the output.
If the 'composite' keyword is set to True, then a `~sunpy.map.CompositeMap` object is returned.
This is useful if the maps are of a different type (e.g. different instruments).
For example, to create a simple Composite Map:

.. code-block:: python

    >>> my_maps = sunpy.map.Map(sunpy.data.sample.EIT_195_IMAGE, sunpy.data.sample.RHESSI_IMAGE, composite=True)  # doctest: +REMOTE_DATA

A `~sunpy.map.CompositeMap` is different from a regular `~sunpy.map.GenericMap` object and therefore different associated methods.
To list which maps are part of your Composite Map use:

.. code-block:: python

    >>> my_maps.list_maps()  # doctest: +REMOTE_DATA
    [<class 'sunpy.map.sources.soho.EITMap'>, <class 'sunpy.map.sources.rhessi.RHESSIMap'>]

The following two examples demonstrate how to create a composite map of AIA and HMI data and how to overlay HMI contours on an AIA map (without creating a composite map object):

* :ref:`sphx_glr_generated_gallery_map_composite_map_AIA_HMI.py`

* :ref:`sphx_glr_generated_gallery_map_hmi_contours_wcsaxes.py`

For a more advanced tutorial on combining data from several maps, see :ref:`sphx_glr_generated_gallery_map_transformations_reprojection_aia_euvi_mosaic.py`.
