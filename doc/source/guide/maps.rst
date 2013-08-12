====
Maps
====

Maps in SunPy are dimensionally-aware data arrays. 
In other words, they are 2-dimensional data associated with a coordinate system. 
In this guide, we will cover some of the basic functionality of maps. 
Once you've read through this guide check out the :doc:`code reference</reference/index>` for a more thorough look at SunPy maps.

------------
Data Support
------------
The map object currently supports the following data sources

- SDO/AIA, SDO/HMI
- STEREO/EUVI, STEREO/COR
- Hinode/XRT
- SOHO/EIT, SOHO/LASCO, SOHO/MDI
- PROBA2/SWAP
- Yohkoh/SXT

1. Creating maps
----------------
SunPy contains a number of example fits files. 
To make things easy, SunPy includes several example files which are used throughout the docs. 
These files have names like `sunpy.AIA_171_IMAGE` and `sunpy.RHESSI_IMAGE`.
To create the sample AIA map type the following into your interactive Python shell::

    import sunpy
    my_map = sunpy.Map(sunpy.AIA_171_IMAGE)

The variable my_map is a SunPy Map object. To create a SunPy Map object from a local fits file try something like the following ::

    my_map = sunpy.Map('/mydirectory/mymap.fits')

SunPy automatically detects the type of file (e.g. fits), what instrument it is 
associated with (e.g. AIA, EIT, LASCO) and will automatically look in the appropriate places for the fits
keywords it needs to interpret the coordinate system. If the type of fits file 
is not recognized then SunPy will try some default fits keywords and return a GenericMap but results
may vary. SunPy can also create maps from the jpg2000 files from
`helioviewer.org <http://helioviewer.org/>`.

2. Creating Custom Maps
-----------------------
It is also possible to create maps using custom data from a simulation for example. To do this you
need to provide Map() with both the data array as well as some basic meta information. If no
header is given then some default values as assumed. Here is a simple example::

    import numpy as np
    data = np.arange(0,100).reshape(10,10)
    header = {'cdelt1': 10, 'cdelt2': 10, 'telescop':'sunpy'}
    my_map = sunpy.Map(data, header)

The format of the header follows the fits standard.

2. Inspecting maps
------------------
A map contains a number of data-associated attributes. To get a quick look at your map simply
type::

    my_map = sunpy.Map(sunpy.AIA_171_IMAGE)
    my_map
    
This will show a representation of the data as well as some of its associated
attributes. A number of other attributes are also available, for example the date, 
exposure time, map center, xrange, yrange
other::

    map_date = my_map.date
    map_exptime = my_map.exposure_time
    map_center = my_map.center
    map_xrange = my_map.xrange
    map_yrange = my_map.yrange
    
To get a list of all of the attributes check the documentation by typing::

    help(my_map)
    
The meta data for the map is accessed by ::

    header = my_map.meta
    
This references the meta data dictionary with the header information as read from the source
file. 

3. Getting at the data
----------------------
The data in a SunPy Map object is accessible through the data attribute. 
Currently, the data is implemented as a NumPy ndarray, so for example, to get 
the 0th element in the array ::

    my_map.data[0,0]
    my_map.data[0][0]
    
One important fact to remember which is intially confusing is that the first index is for the 
y direction while the second index is for the x direction! For more information about indexing 
please refer to the `Numpy documentation <http://www.scipy.org/Tentative_NumPy_Tutorial#head-864862d3f2bb4c32f04260fac61eb4ef34788c4c>`.
Common ndarray attributes, such as shape and dtype, are accessible through the SunPy Map object ::

    my_map.shape
    my_map.dtype

If you'd like to use the data in a SunPy Map object elsewhere, you can use ::

    var = my_map.data
    # or
    var = my_map.data.copy()
    
Basic statistical functions on the data array are also passed through to Map objects::

    my_map.min()
    my_map.max()
    my_map.mean()

4. Creating a plot of your map
------------------------------
The SunPy map object has its own built-in plot methods so that it is easy to
quickly view your map on the screen. To create a plot just type::

    my_map.peek()
    
This will open a matplotlib plot right on your screen.
In addition, to enable users to modify the plot it is possible to grab the
matplotlib figure object by using the plot() command instead of the show() 
command. This makes it possible to use the SunPy plot as the foundation for a 
more complicated figure.

5. Overlaying Maps
------------------
The Map() method described above can also handle a list of maps. If a list in inputs
is supplied, Map() will return a list of maps as the output.  However, if the
'composite' keyword is set to True, then a CompositeMap object is returned.  This is useful if the maps are
of a different type (e.g. different instruments).  For example, to create a simple composite map::

    my_maps = sunpy.Map(sunpy.EIT_195_IMAGE, sunpy.RHESSI_IMAGE, composite=True)

A CompositeMap is different from a regular SunPy Map objectand therefore different associated methods.
To list which maps are part of your composite map use::

    my_maps.list_maps()

Similar to all SunPy data objects, the composite map also has an associated show() method and a 
number of associated methods to customize your plot. For example, the following code turns 
adds a new map (which must be instantiated first), sets its transparency to 25%, turns on contours from 50% to 90% for the second map, 
and then plots the result::

    my_maps.add_map(sunpy.Map(sunpy.AIA_171_IMAGE))
    my_maps.set_alpha(2,0.5)
    my_maps.set_levels(1,[50,60,70,80,90], percent = True)
    my_maps.peek()

This is not a particularly pretty plot but it shows what SunPy can do!

5. Working with your map
------------------------
Part of the philosophy of the map object is to provide most of the basic
functionality that a scientist would want therefore a map also contains a number
of map-specific methods such as resizing a map or grabbing a subview. To get 
a list of the methods available for a map type::

    help(my_map)
    
and check out the methods section!

