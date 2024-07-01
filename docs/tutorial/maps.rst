.. _sunpy-tutorial-maps:

****
Maps
****

In this section of the tutorial, you will learn about the `Map <sunpy.map.GenericMap>` object.
Map objects hold two-dimensional data along with the accompanying metadata.
They can be used with any two-dimensional data array with two spatial axes and FITS-compliant metadata.
As you will see in this tutorial, the primary advantage of using a Map object over a bare NumPy array is the ability to perform coordinate aware operations on the underlying array, such as rotating the Map to remove the roll angle of an instrument or cropping a Map to a specific field of view.
Additionally, because Map is capable of ingesting data from many different sources, it provides a common interface for any two-dimensional data product.

By the end of this tutorial, you will learn how to create a Map, access the underlying data and metadata, convert between physical coordinates and pixel coordinates of a Map, and the basics of visualizing a Map.
Additionally, you will learn how to combine multiple Maps into a `~sunpy.map.MapSequence` or a `~sunpy.map.CompositeMap`.

.. note::

    In this section and in :ref:`sunpy-tutorial-timeseries`, we will use the sample data included with sunpy.
    These data are primarily useful for demonstration purposes or simple debugging.
    These files have names like ``sunpy.data.sample.AIA_171_IMAGE`` and ``sunpy.data.sample.RHESSI_IMAGE`` and are automatically downloaded to your computer as you need them.
    Once downloaded, these sample data files will be paths to their location on your computer.

.. _sunpy-tutorial-map-creating-maps:

Creating Maps
=============

To create a `~sunpy.map.Map` from a sample Atmospheric Imaging Assembly (AIA) image,

.. code-block:: python

    >>> import sunpy.map
    >>> import sunpy.data.sample  # doctest: +REMOTE_DATA

    >>> sunpy.data.sample.AIA_171_IMAGE  # doctest: +REMOTE_DATA
    PosixPath('.../AIA20110607_063302_0171_lowres.fits')
    >>> my_map = sunpy.map.Map(sunpy.data.sample.AIA_171_IMAGE)  # doctest: +REMOTE_DATA

In many cases sunpy automatically detects the type of file as well as the the instrument associated with it.
In this case, we have a FITS file from the AIA instrument onboard the Solar Dynamics Observatory (SDO).
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

.. _sunpy-tutorial-map-inspecting-maps:

Inspecting Map Metadata
=======================

The metadata for a Map is exposed via attributes on the Map.
These attributes can be accessed by typing ``my_map.<attribute-name>``.
For example, to access the date of the observation,

.. code-block:: python

    >>> my_map.date  # doctest: +REMOTE_DATA
    <Time object: scale='utc' format='isot' value=2011-06-07T06:33:02.770>

Notice that this is an `~astropy.time.Time` object which we discussed in the previous :ref:`sunpy-tutorial-times` section of the tutorial.
Similarly, we can access the exposure time of the image,

.. code-block:: python

    >>> my_map.exposure_time  # doctest: +REMOTE_DATA
    <Quantity 0.234256 s>

Notice that this returns an `~astropy.units.Quantity` object which we discussed in the previous :ref:`sunpy-tutorial-units` section of the tutorial.
The full list of attributes can be found in the reference documentation for `~sunpy.map.GenericMap`.
These metadata attributes are all derived from the underlying FITS metadata, but are represented as rich Python objects, rather than simple strings or numbers.

.. _sunpy-tutorial-map-data:

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
    np.float32(-95.92475)

The first index corresponds to the y direction and the second to the x direction in the two-dimensional pixel coordinate system.
For more information about indexing, please refer to the `numpy documentation <https://numpy.org/doc/stable/user/basics.indexing.html#indexing-on-ndarrays>`__.

Data attributes like dimensionality and type are also accessible as attributes on ``my_map``:

.. code-block:: python

    >>> my_map.dimensions  # doctest: +REMOTE_DATA
    PixelPair(x=<Quantity 1024. pix>, y=<Quantity 1024. pix>)
    >>> my_map.dtype  # doctest: +REMOTE_DATA
    dtype('float32')

Additionally, there are several methods that provide basic summary statistics of the data:

.. code-block:: python

    >>> my_map.min()  # doctest: +REMOTE_DATA
    np.float32(-129.78036)
    >>> my_map.max()  # doctest: +REMOTE_DATA
    np.float32(192130.17)
    >>> my_map.mean()  # doctest: +REMOTE_DATA
    np.float32(427.02252)

.. _sunpy-tutorial-map-coordinates-wcs:

Coordinates, and the World Coordinate System
============================================

In :ref:`sunpy-tutorial-coordinates`, you learned how to define coordinates with `~astropy.coordinates.SkyCoord` using different solar coordinate frames.
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

This tells us the location of the spacecraft, in this case SDO, when it recorded this particular observation, as derived from the FITS metadata.

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
Each Map has an associated World Coordinate System, or WCS, which is derived from the underlying metadata and expressed as an `astropy.wcs.WCS` object.
The WCS is accessible as an attribute:

.. code-block:: python

    >>> my_map.wcs  # doctest: +REMOTE_DATA
    WCS Keywords
    <BLANKLINE>
    Number of WCS axes: 2
    CTYPE : 'HPLN-TAN' 'HPLT-TAN'
    CRVAL : np.float64(0.00089530541880571) np.float64(0.00038493926472939)
    CRPIX : np.float64(512.5) np.float64(512.5)
    PC1_1 PC1_2  : np.float64(0.99999706448085) np.float64(0.0024230207763071)
    PC2_1 PC2_2  : np.float64(-0.0024230207763071) np.float64(0.99999706448085)
    CDELT : np.float64(0.00066744222222222) np.float64(0.00066744222222222)
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

.. _sunpy-tutorial-map-plotting-maps:

Visualizing Maps
================

.. plot::
    :nofigs:
    :context: close-figs
    :show-source-link: False

    # This is here to put my_map in the scope of the plot directives.
    # This avoids repeating code in the example source code that is actually displayed.
    # This snippet of code is not visible in the rendered documentation.
    import sunpy.map
    import sunpy.data.sample
    from astropy.coordinates import SkyCoord
    import astropy.units as u
    my_map = sunpy.map.Map(sunpy.data.sample.AIA_171_IMAGE)

In the :ref:`sunpy-tutorial-map-creating-maps` section, you learned how to generate a quicklook summary of a Map.
However, the Map object also has a :meth:`~sunpy.map.GenericMap.plot` method that allows for more fine-grained control over how the Map is visualized and is especially useful for generating publication-quality plots.
In this section of the tutorial, you will learn how to build up an increasingly detailed visualization of a Map, including adjusting the colormap and normalization and and overlaying coordinates and contours.

Basic Plotting
--------------

First, let's create a basic plot of our Map, including a colorbar,

.. plot::
    :include-source:
    :context: close-figs

    import matplotlib.pyplot as plt

    fig = plt.figure()
    ax = fig.add_subplot(projection=my_map)
    my_map.plot(axes=ax)
    plt.colorbar()
    plt.show()

.. note::

    We imported `matplotlib.pyplot` in order to create the figure and the axis we plotted on our map onto.
    Under the hood, sunpy uses `matplotlib` to visualize the image meaning that plots built with sunpy can be further customized using `matplotlib`.
    **However, for the purposes of this tutorial, you do not need to be familiar with Matplotlib.**
    For a series of detailed examples showing how to customize your Map plots, see the :ref:`Plotting section of the Example Gallery <sphx_glr_generated_gallery_plotting>` as well as the documentation for `astropy.visualization.wcsaxes`.

Note that the title and colormap have been set by sunpy based on the observing instrument and wavelength.
Furthermore, the tick and axes labels have been automatically set based on the coordinate system of the Map.

Looking at the plot above, you likely notice that the resulting image is a bit dim.
To fix this, we can use the ``clip_interval`` keyword to automatically adjust the colorbar limits to clip out the dimmest 1% and the brightest 0.5% of pixels.

.. plot::
    :include-source:
    :context: close-figs

    fig = plt.figure()
    ax = fig.add_subplot(projection=my_map)
    my_map.plot(axes=ax, clip_interval=(1, 99.5)*u.percent)
    plt.colorbar()
    plt.show()

Changing the Colormap and Normalization
---------------------------------------

Historically, particular colormaps are assigned to images based on what instrument they are from and what wavelength is being observed.
By default, sunpy will select the colormap based on the available metadata.
This default colormap is available as an attribute,

.. code-block:: python

    >>> my_map.cmap.name  # doctest: +REMOTE_DATA
    'sdoaia171'

When visualizing a Map, you can change the colormap using the ``cmap`` keyword argument.
For example, you can use the 'inferno' colormap from `matplotlib`:

.. plot::
    :include-source:
    :context: close-figs

    fig = plt.figure()
    ax = fig.add_subplot(projection=my_map)
    my_map.plot(axes=ax, cmap='inferno', clip_interval=(1,99.5)*u.percent)
    plt.colorbar()
    plt.show()

.. note::

    sunpy provides specific colormaps for many different instruments.
    For a list of all colormaps provided by sunpy, see the documentation for `sunpy.visualization.colormaps`.

The normalization, or the mapping between the data values and the colors in our colormap, is also determined based on the underlying metadata.
Notice that in the plots we've made so far, the ticks on our colorbar are not linearly spaced.
Just like in the case of the colormap, we can use a normalization other than the default by passing a keyword argument to the :meth:`~sunpy.map.GenericMap.plot` method.
For example, we can use a logarithmic normalization instead:

.. plot::
    :include-source:
    :context: close-figs

    import matplotlib.colors

    fig = plt.figure()
    ax = fig.add_subplot(projection=my_map)
    my_map.plot(norm=matplotlib.colors.LogNorm())
    plt.colorbar()
    plt.show()

.. note::

    You can also view or make changes to the default settings through the ``sunpy.map.GenericMap.plot_settings`` dictionary.
    See :ref:`sphx_glr_generated_gallery_plotting_map_editcolormap.py` for an example of of how to change the default plot settings.

.. _sunpy-tutorial-map-wcsaxes-plotting:

Overlaying Contours and Coordinates
-----------------------------------

When plotting images, we often want to highlight certain features or overlay certain data points.
There are several methods attached to Map that make this task easy.
For example, we can draw contours around the brightest 0.5% percent of pixels in the image:

.. plot::
    :include-source:
    :context: close-figs

    fig = plt.figure()
    ax = fig.add_subplot(projection=my_map)
    my_map.plot(axes=ax, clip_interval=(1,99.5)*u.percent)
    my_map.draw_contours([2, 5, 10, 50, 90] * u.percent, axes=ax)
    plt.show()

Additionally, the solar limb, as determined by the location of the observing instrument at the time of the observation, can be easily overlaid on an image:

.. plot::
    :include-source:
    :context: close-figs

    fig = plt.figure()
    ax = fig.add_subplot(projection=my_map)
    my_map.plot(axes=ax, clip_interval=(1,99.5)*u.percent)
    my_map.draw_limb(axes=ax, color='C0')
    plt.show()

We can also overlay a box denoting a particular a region of interest as expressed in world coordinates using the the coordinate frame of our image:

.. plot::
    :include-source:
    :context: close-figs

    roi_bottom_left = SkyCoord(Tx=-300*u.arcsec, Ty=-100*u.arcsec, frame=my_map.coordinate_frame)
    roi_top_right = SkyCoord(Tx=200*u.arcsec, Ty=400*u.arcsec, frame=my_map.coordinate_frame)
    fig = plt.figure()
    ax = fig.add_subplot(projection=my_map)
    my_map.plot(axes=ax, clip_interval=(1,99.5)*u.percent)
    my_map.draw_quadrangle(roi_bottom_left, top_right=roi_top_right, axes=ax, color='C0')
    plt.show()

Because our visualization knows about the coordinate system of our Map, it can transform any coordinate to the coordinate frame of our Map and then use the underlying WCS that we discussed in the :ref:`sunpy-tutorial-map-coordinates-wcs` section to translate this to a pixel position.
This makes it simple to plot *any* coordinate on top of our Map using the :meth:`~astropy.visualization.wcsaxes.WCSAxes.plot_coord` method.
The following example shows how to plot some points on our Map, including the center coordinate of our Map:

.. plot::
    :include-source:
    :context: close-figs

    coords = SkyCoord(Tx=[100,1000] * u.arcsec, Ty=[100,1000] * u.arcsec, frame=my_map.coordinate_frame)

    fig = plt.figure()
    ax = fig.add_subplot(projection=my_map)
    my_map.plot(axes=ax, clip_interval=(1,99.5)*u.percent)
    ax.plot_coord(coords, 'o')
    ax.plot_coord(my_map.center, 'X')
    plt.show()

.. note::

    Map visualizations can be heavily customized using both `matplotlib` and `astropy.visualization.wcsaxes`.
    See the :ref:`Plotting section of the Example Gallery <sphx_glr_generated_gallery_plotting>` for more detailed examples of how to customize Map visualizations.

.. _sunpy-tutorial-map-cropping-maps:

Cropping Maps and Combining Pixels
==================================

In analyzing images of the Sun, we often want to choose a smaller portion of the full disk to look at more closely.
Let's use the region of interest we defined above to crop out that portion of our image:

.. plot::
    :include-source:
    :context: close-figs

    my_submap = my_map.submap(roi_bottom_left, top_right=roi_top_right)

    fig = plt.figure()
    ax = fig.add_subplot(projection=my_submap)
    my_submap.plot(axes=ax)
    plt.show()

Additionally, we also may want to combine multiple pixels into a single pixel (a "superpixel") to, for example, increase our signal-to-noise ratio.
We can accomplish this with the `~sunpy.map.GenericMap.superpixel` method by specifying how many pixels, in each dimension, we want our new superpixels to contain.
For example, we can combine 4 pixels in each dimension such that our new superpixels contain 16 original pixels:

.. plot::
    :include-source:
    :context: close-figs

    my_super_submap = my_submap.superpixel((5,5)*u.pixel)

    fig = plt.figure()
    ax = fig.add_subplot(projection=my_super_submap)
    my_super_submap.plot(axes=ax)
    plt.show()

.. note::

    Map provides additional methods for manipulating the underlying image data.
    See the reference documentation for `~sunpy.map.GenericMap` for a complete list of available methods as well as the :ref:`Map section of the Example Gallery <sphx_glr_generated_gallery_map>` for more detailed examples.

.. _sunpy-tutorial-map-map-sequences:

Map Sequences
=============

While `~sunpy.map.GenericMap` can only contain a two-dimensional array and metadata corresponding to a single observation, a `~sunpy.map.MapSequence` is comprised of an ordered list of maps.
By default, the Maps are ordered by their observation date, from earliest to latest date.
A `~sunpy.map.MapSequence` can be created by supplying multiple existing maps:

.. code-block:: python

    >>> another_map = sunpy.map.Map(sunpy.data.sample.EIT_195_IMAGE)  # doctest: +REMOTE_DATA
    >>> map_seq = sunpy.map.Map([my_map, another_map], sequence=True)  # doctest: +REMOTE_DATA

A map sequence can be indexed in the same manner as a list.
For example, the following returns the same information as in :ref:`sunpy-tutorial-map-creating-maps`:

.. code-block:: python

    >>> map_seq.maps[0]   # doctest: +REMOTE_DATA
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

MapSequences can hold maps that have different shapes.
To test if all the maps in a `~sunpy.map.MapSequence` have the same shape:

.. code-block:: python

    >>> map_seq.all_maps_same_shape()  # doctest: +REMOTE_DATA
    np.True_

It is often useful to return the image data in a `~sunpy.map.MapSequence` as a single three dimensional NumPy `~numpy.ndarray`:

.. code-block:: python

    >>> map_seq_array = map_seq.as_array()  # doctest: +REMOTE_DATA

Since all of the maps in our sequence of the same shape, the first two dimensions of our combined array will be the same as the component maps while the last dimension will correspond to the number of maps in the map sequence.
We can confirm this by looking at the shape of the above array.

.. code-block:: python

    >>> map_seq_array.shape  # doctest: +REMOTE_DATA
    (1024, 1024, 2)

.. warning::

    `~sunpy.map.MapSequence` does not automatically perform any coalignment between the maps comprising a sequence.
    For information on coaligning images and compensating for solar rotation, see :ref:`this section of the Example Gallery <sphx_glr_generated_gallery_map_transformations>` as well as the `sunkit_image.coalignment` module.
