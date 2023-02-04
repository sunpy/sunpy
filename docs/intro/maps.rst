.. _map_guide:

****
Maps
****

Map objects in **sunpy** are two-dimensional data associated with a coordinate system.
This class offers many advantages over using packages such as `astropy.io.fits` to open 2D solar data in Python.
**sunpy** Maps be used with any image with two spatial axes that has FITS standard compliant metadata.
There are also 19 sources (see :ref:`map` for a complete list) that provide instrument-specific information when plotting and in some cases interpret missing metadata when it is missing from the original file.

Part of the philosophy of the Map object is to provide most of the basic functionality that a scientist would want.
Therefore, a Map also contains a number of methods such as resizing or grabbing a subview.
See `~sunpy.map.mapbase.GenericMap` for summaries of all attributes and methods.

:ref:`creating-maps` and :ref:`inspecting-maps` demonstrate how to create a Map object and view a summary of the data in the object.
:ref:`map-data` describes how data is stored and accessed in a Map, and :ref:`plotting-maps` introduces some basics of visualizing Map data in **sunpy**.
:ref:`map-sequences` and :ref:`composite-maps` detail how **sunpy** can handle multiple Maps from similar and different instruments.
Lastly, :ref:`custom-maps` shows how you can create Maps from data sources not supported in **sunpy**.

.. _creating-maps:

Creating Maps
=============

To make things easy, **sunpy** can download sample data which are used throughout the docs.
These files have names like ``sunpy.data.sample.AIA_171_IMAGE`` and ``sunpy.data.sample.RHESSI_IMAGE``.
To create a `~sunpy.map.Map` from a sample AIA image, type the following into your Python shell:

.. code-block:: python

    >>> import sunpy
    >>> import sunpy.map
    >>> import sunpy.data.sample  # doctest: +REMOTE_DATA

    >>> my_map = sunpy.map.Map(sunpy.data.sample.AIA_171_IMAGE)  # doctest: +REMOTE_DATA

The variable ``my_map`` is an `~sunpy.map.sources.AIAMap` object.
To create a Map from a local FITS file try the following:

.. code-block:: python

    >>> my_map = sunpy.map.Map('/mydirectory/mymap.fits')   # doctest: +SKIP

**sunpy** automatically detects the type of file (e.g. FITS), the instrument associated with it (e.g. AIA, EIT, LASCO), and finds the FITS keywords it needs to interpret the coordinate system.
If the type of FITS file is not recognized then **sunpy** will try some default FITS keywords and return a `~sunpy.map.GenericMap`, but results may vary.
**sunpy** can also create Maps using the JPEG 2000 files from `helioviewer.org <https://helioviewer.org/>`__.

The :ref:`sphx_glr_generated_gallery_acquiring_data` section of the Example Gallery has several examples of retrieving data using the `~sunpy.net.Fido` tool.
For more details, check out the :ref:`fido_guide` guide.
See :ref:`sphx_glr_generated_gallery_saving_and_loading_data_genericmap_in_fits.py` for a demonstration of saving a map as a FITS file and `~sunpy.map.GenericMap.save` for more details on saving.

.. _inspecting-maps:

Inspecting Maps
===============

A Map contains a number of data-associated attributes.
To get a quick look at your Map simply type:

.. code-block:: python

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

This will show a representation of the data as well as some of its associated attributes.
Typing the above command in a Jupyter Notebook will show a rich HTML view of the table along with two plots of your data.
The HTML view can also be accessed using the :func:`~sunpy.map.GenericMap.quicklook` method, which will open the view in your default browser.

A number of other attributes are also available.
For example, the `~sunpy.map.GenericMap.date`, `~sunpy.map.GenericMap.exposure_time`, `~sunpy.map.GenericMap.center` and others (see `~sunpy.map.GenericMap`).
The full list can be found on `~sunpy.map.GenericMap`:

.. code-block:: python

    >>> map_date = my_map.date  # doctest: +REMOTE_DATA
    >>> map_exptime = my_map.exposure_time  # doctest: +REMOTE_DATA
    >>> map_center = my_map.center  # doctest: +REMOTE_DATA

To get a list of all of the attributes check the documentation by typing:

.. code-block:: python

    >>> help(my_map)  # doctest: +SKIP

Many attributes and functions of the map classes accept and return `~astropy.units.quantity.Quantity` or `~astropy.coordinates.SkyCoord` objects.
Please refer to :ref:`units-sunpy` and :ref:`coordinates-sunpy` for more details.

The metadata for the map is accessed by:

.. code-block:: python

    >>> header = my_map.meta  # doctest: +REMOTE_DATA

This references the metadata dictionary with the header information as read from the source file.
To see if the metadata of a Map source has been modified, see :ref:`sphx_glr_generated_gallery_map_map_metadata_modification.py` for a demonstration.

.. _map-data:

Map Data
========

The data in a Map object is accessible through the `~sunpy.map.GenericMap.data` attribute.
The data is stored as a NumPy `~numpy.ndarray`.
For example, to get the 0th element in the array:

.. code-block:: python

    >>> my_map.data[0, 0]  # doctest: +REMOTE_DATA
    -95.92475
    >>> my_map.data[0][0]  # doctest: +REMOTE_DATA
    -95.92475

The first index is for the y direction while the second index is for the x direction.
For more information about indexing, please refer to the `numpy documentation <https://numpy.org/doc/stable/user/basics.indexing.html#indexing-on-ndarrays>`__.

Data attributes like `~numpy.ndarray.dtype` and `~sunpy.map.GenericMap.dimensions` are accessible through a GenericMap object:

.. code-block:: python

    >>> my_map.dimensions  # doctest: +REMOTE_DATA
    PixelPair(x=<Quantity 1024. pix>, y=<Quantity 1024. pix>)
    >>> my_map.dtype  # doctest: +REMOTE_DATA
    dtype('float32')

Here, the dimensions attribute is similar to the `~numpy.ndarray.shape` attribute, however returning an `~astropy.units.quantity.Quantity`.

You can store the data of a `~sunpy.map.GenericMap` object in a separate `~numpy.ndarray` by either of the following actions:

.. code-block:: python

    >>> var = my_map.data  # doctest: +REMOTE_DATA
    >>> var = my_map.data.copy()  # doctest: +REMOTE_DATA

To create a complete copy of a Map object that is entirely independent of the original, use the built-in `copy.deepcopy` function:

.. code-block:: python

    >>> import copy   # doctest: +REMOTE_DATA
    >>> my_map_deepcopy = copy.deepcopy(my_map)   # doctest: +REMOTE_DATA

A deepcopy ensures that any changes in the original Map object are not reflected in the copied object and vice versa.
Note that this copies the data of the Map object as well as all of the other attributes and methods.

Some basic statistical functions are built into Map objects:

.. code-block:: python

    >>> my_map.min()  # doctest: +REMOTE_DATA
    -129.78036
    >>> my_map.max()  # doctest: +REMOTE_DATA
    192130.17
    >>> my_map.mean()  # doctest: +REMOTE_DATA
    427.02252

All the other `~numpy.ndarray` functions and attributes can be accessed through the data array directly.
For example:

.. code-block:: python

    >>> my_map.data.std()  # doctest: +REMOTE_DATA
    826.41016

.. _plotting-maps:

Plotting Maps
=============

The `~sunpy.map.GenericMap` object has a built-in plot method such that it is easy to quickly view your map.
To create a plot just type:

.. code-block:: python

    >>> my_map.peek()   # doctest: +SKIP

This will open a Matplotlib plot on your screen.
In addition, it is possible to grab the Matplotlib Axes object by using the `~sunpy.map.GenericMap.plot()` command.
This makes it possible to use the **sunpy** plot as the foundation for a more complicated figure.
For more information about this and some examples see :ref:`plotting`.
Check out the following foundational examples in the Example Gallery for plotting Maps:

* :ref:`sphx_glr_generated_gallery_plotting_aia_example.py`

* :ref:`sphx_glr_generated_gallery_plotting_wcsaxes_plotting_example.py`

* :ref:`sphx_glr_generated_gallery_plotting_map_editcolormap.py`

* :ref:`sphx_glr_generated_gallery_plotting_grid_plotting.py`

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
    smap = sunpy.map.Map(sunpy.data.sample.AIA_171_IMAGE)
    fig = plt.figure()
    smap.plot(cmap=plt.cm.Greys_r)
    plt.colorbar()
    plt.show()

You can also view or make changes to the default settings through the ``sunpy.map.GenericMap.plot_settings`` dictionary.
See :ref:`sphx_glr_generated_gallery_plotting_map_editcolormap.py` for an example of this workflow for changing plot settings.


Colormaps and Normalization
---------------------------

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

    smap = sunpy.map.Map(sunpy.data.sample.AIA_171_IMAGE)
    cmap = plt.get_cmap('sohoeit171')

    fig = plt.figure()
    smap.plot(cmap=cmap)
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

    smap = sunpy.map.Map(sunpy.data.sample.AIA_171_IMAGE)

    fig = plt.figure(figsize=(4, 9))

    ax1 = fig.add_subplot(2, 1, 1, projection=smap)
    smap.plot(norm=colors.Normalize(), title='Linear normalization')
    plt.colorbar()

    ax2 = fig.add_subplot(2, 1, 2, projection=smap)
    smap.plot(norm=colors.LogNorm(), title='Logarithmic normalization')
    plt.colorbar()

    plt.show()

Note how the colorbar does not change since these two plots share the same colormap.
Meanwhile, the data values associated with each color do change because the normalization is different.


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

or by providing a directory full of image files:

.. code-block:: python

    >>> mc = sunpy.map.Map('path/to/my/files/*.fits', sequence=True)   #  doctest: +SKIP

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

.. _custom-maps:

Creating Custom Maps
====================

It is also possible to create Maps using custom data (e.g. from a simulation or an observation from a data source that is not explicitly supported in **sunpy**).
To do this, you need to provide `sunpy.map.Map` with both the data array as well as appropriate meta information.
The meta information informs `sunpy.map.Map` of the correct coordinate information associated with the data array and should be provided to `sunpy.map.Map` in the form of a header as a `dict` or `~sunpy.util.MetaDict`.
See this :ref:`sphx_glr_generated_gallery_map_map_from_numpy_array.py` for a brief demonstration of generating a Map from a data array.

The keys required for the header information follow the `FITS standard <https://fits.gsfc.nasa.gov/fits_dictionary.html>`__.
All the meta keywords used when creating a Map along with their description are listed below.

`cunit1`: Units of the coordinate increments along naxis1 e.g. arcsec \**required
`cunit2`: Units of the coordinate increments along naxis2 e.g. arcsec \**required
`crval1`: Coordinate value at reference point on naxis1 \**required
`crval2`: Coordinate value at reference point on naxis2 \**required
`cdelt1`: Spatial scale of pixels for naxis1, i.e. coordinate increment at reference point
`cdelt2`: Spatial scale of pixels for naxis2, i.e. coordinate increment at reference point
`crpix1`: Pixel coordinate at reference point naxis1
`crpix2`: Pixel coordinate at reference point naxis2
`ctype1`: Coordinate type projection along naxis1 of data e.g. HPLT-TAN
`ctype2`: Coordinate type projection along naxis2 of data e.g. HPLN-TAN
`hgln_obs`: Heliographic longitude of observation
`hglt_obs`: Heliographic latitude of observation
`dsun_obs`: distance to Sun from observation in metres
`rsun_obs`: radius of Sun in meters from observation
`dateobs`: date of observation e.g. 2013-10-28 00:00
`date_obs`: date of observation e.g. 2013-10-28 00:00
`rsun_ref`: reference radius of Sun in meters
`solar_r`: radius of Sun in meters from observation
`radius`: radius of Sun in meters from observation
`crln_obs`: Carrington longitude of observation
`crlt_obs`: Heliographic latitude of observation
`solar_b0`: Solar B0 angle
`detector`: name of detector e.g. AIA
`exptime`: exposure time of observation, in seconds e.g 2
`instrume`: name of instrument
`wavelnth`: wavelength of observation
`waveunit`: unit for which observation is taken e.g. angstom
`obsrvtry`: name of observatory of observation
`telescop`: name of telescope of observation
`lvl_num`: FITS processing level
`crota2`: Rotation of the horizontal and vertical axes in degrees
`PC1_1`: Matrix element PCi_j describing the rotation required to align solar North with the top of the image.
`PC1_2`: Matrix element PCi_j describing the rotation required to align solar North with the top of the image.
`PC2_1`: Matrix element PCi_j describing the rotation required to align solar North with the top of the image.
`PC2_2`: Matrix element PCi_j describing the rotation required to align solar North with the top of the image.
`CD1_1`: Matrix element CDi_j describing the rotation required to align solar North with the top of the image.
`CD1_2`: Matrix element CDi_j describing the rotation required to align solar North with the top of the image.
`CD2_1`: Matrix element CDi_j describing the rotation required to align solar North with the top of the image.
`CD2_2`: Matrix element CDi_j describing the rotation required to align solar North with the top of the image.

``sunpy`` provides a helper, :func:`~sunpy.map.header_helper.make_fitswcs_header`, to assist in creating a header that contains the correct metadata.
This will return a header with the appropriate FITS keywords once the Map data array and an `astropy.coordinates.SkyCoord` or `sunpy.coordinates.frames` is provided.
The `astropy.coordinates.SkyCoord` is defined by the user and contains information on the reference frame, reference coordinate, and observer location.
This function returns a `sunpy.util.MetaDict`.
The `astropy.coordinates.SkyCoord` or `sunpy.coordinates.frames` must contain an observation time.

The :func:`~sunpy.map.header_helper.make_fitswcs_header` function also takes optional keyword arguments including ``reference_pixel`` and ``scale`` that describe the pixel coordinate at the reference coordinate (defined by the `~astropy.coordinates.SkyCoord`) and the spatial scale of the pixels, respectively.
If neither of these are given their values default to the center of the data array and 1 arcsec, respectively.

Here's an example of creating a header from some generic data and an `astropy.coordinates.SkyCoord`:

.. code-block:: python

    >>> import numpy as np
    >>> from astropy.coordinates import SkyCoord
    >>> import astropy.units as u

    >>> from sunpy.coordinates import frames
    >>> from sunpy.map.header_helper import make_fitswcs_header

    >>> data = np.arange(0,100).reshape(10,10)
    >>> coord = SkyCoord(0*u.arcsec, 0*u.arcsec, obstime = '2013-10-28', observer = 'earth', frame = frames.Helioprojective)
    >>> header = make_fitswcs_header(data, coord)
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
    lonpole: 180.0
    latpole: 0.0
    mjdref: 0.0
    date-obs: 2013-10-28T00:00:00.000
    rsun_ref: 695700000.0
    dsun_obs: 148644585949.49
    hgln_obs: 0.0
    hglt_obs: 4.7711570596394
    naxis: 2
    naxis1: 10
    naxis2: 10
    pc1_1: 1.0
    pc1_2: -0.0
    pc2_1: 0.0
    pc2_2: 1.0
    rsun_obs: 965.3829548285768

From this we can see now that the function returned a `sunpy.util.MetaDict` that populated the standard FITS keywords with information provided by the passed `astropy.coordinates.SkyCoord`, and the data array.
Since the ``reference_pixel`` and keywords were not passed in the example above, the values of ``crpix`` and ``cdelt`` were set to the default values.

These keywords can be passed to the function in the form of an `astropy.units.Quantity` with associated units.
Here's another example of passing ``reference_pixel`` and ``scale`` to the function:

.. code-block:: python

    >>> header = make_fitswcs_header(data, coord,
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
    lonpole: 180.0
    latpole: 0.0
    mjdref: 0.0
    date-obs: 2013-10-28T00:00:00.000
    rsun_ref: 695700000.0
    dsun_obs: 148644585949.49
    hgln_obs: 0.0
    hglt_obs: 4.7711570596394
    naxis: 2
    naxis1: 10
    naxis2: 10
    pc1_1: 1.0
    pc1_2: -0.0
    pc2_1: 0.0
    pc2_2: 1.0
    rsun_obs: 965.3829548285768

As we can see, a list of WCS and observer meta information is contained within the generated headers, however we may want to include other meta information including the observatory name, the wavelength and waveunit of the observation.
Any of the keywords in the dictionary returned by :func:`~sunpy.map.header_helper.meta_keywords` can be passed to the :func:`~sunpy.map.header_helper.make_fitswcs_header` and will then populate the returned MetaDict header.
Furthermore, the following observation keywords can be passed to the `~sunpy.map.header_helper.make_fitswcs_header` function: ``observatory``, ``instrument``, ``telescope``, ``wavelength``, ``exposure``.

An example of creating a header with these additional keywords:

.. code-block:: python

    >>> header = make_fitswcs_header(data, coord,
    ...                                        reference_pixel = u.Quantity([5, 5]*u.pixel),
    ...                                        scale = u.Quantity([2, 2] *u.arcsec/u.pixel),
    ...                                        telescope = 'Test case', instrument = 'UV detector',
    ...                                        wavelength = 1000*u.angstrom)
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
    lonpole: 180.0
    latpole: 0.0
    mjdref: 0.0
    date-obs: 2013-10-28T00:00:00.000
    rsun_ref: 695700000.0
    dsun_obs: 148644585949.49
    hgln_obs: 0.0
    hglt_obs: 4.7711570596394
    instrume: UV detector
    telescop: Test case
    wavelnth: 1000.0
    waveunit: Angstrom
    naxis: 2
    naxis1: 10
    naxis2: 10
    pc1_1: 1.0
    pc1_2: -0.0
    pc2_1: 0.0
    pc2_2: 1.0
    rsun_obs: 965.3829548285768

From these header MetaDict's that are generated, we can now create a custom map:

.. code-block:: python

    >>> my_map = sunpy.map.Map(data, header)
