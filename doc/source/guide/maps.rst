====
Maps
====

Maps in SunPy are spatially-aware data arrays. In other words, they are 2-dimensional data associated with a coordinate system. In this guide, we will 
cover some of the basic functionality of maps. Once you've read through this guide check out the :doc:`code reference</reference/index>` for a more
thorough look at SunPy maps.

------------
Data Support
------------
The map object currently supports the following data sources

- SDO/AIA
- STEREO/EUVI, STEREO/COR
- Hinode/XRT
- SOHO/EIT, SOHO/LASCO, SOHO/MDI
- PROBA2/SWAP
- Yohkoh/SXT

1. Creating maps
----------------
SunPy contains a number of example fits files. To make things easy,
SunPy includes several example files which are used throughout the docs. These
files have names like `sunpy.AIA_171_IMAGE` and `sunpy.RHESSI_IMAGE`.
To create the sample AIA map type the following into your interactive Python shell::

    import sunpy
    my_map = sunpy.make_map(sunpy.AIA_171_IMAGE)

The variable my_map is a SunPy map. To create a map from a local fits file try
something like the following ::

    my_map = sunpy.make_map('/mydirectory/mymap.fits')

SunPy should automatically detect the type of file (e.g. fits), what instrument it is 
associated with (e.g. AIA, EIT, LASCO) and will automatically look in the appropriate places for the fits
keywords it needs to interpret the coordinate system. If the type of fits file 
is not recognized then SunPy will try some default fits keywords but results
may vary. The list of files which are currently supported by SunPy can be found by looking at the 
documentation for make_map(). SunPy can also create maps from the jpg2000 files from
`helioviewer.org <http://helioviewer.org/>`.

2. Creating Custom Maps
-----------------------
It is also possible to create maps using custom data from a simulation for example. To do this you
need to provide make_map() with both the data array as well as some basic header information. If no
header is given then some default values as assumed. Here is a simple example::

    import numpy as np
    data = np.arange(0,100).reshape(10,10)
    header = {'cdelt1': 10, 'cdelt2': 10, 'telescop':'sunpy'}
    my_map = sunpy.make_map(data, header)

The format of the header follows the fits standard.

2. Inspecting maps
------------------
A map contains a number of data-associated attributes. To get a quick look at your map simply
type::

    my_map = sunpy.make_map(sunpy.AIA_171_IMAGE)
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
    
The header of the map can be accessed by using the following method::

    header = my_map.get_header()
    
This returns a dictionary with the header information as read from the source
file. 

3. Getting at the data
----------------------
As a SunPy map inherits from the NumPy ndarray it can be accessed in the same
way, through indexing. For example, to get the 0th element in the array ::

    my_map[0,0]
    my_map[0][0]
    
One important fact to remember which is intially confusing is that the first index is for the 
y direction while the second index is for the x direction! For more information about indexing 
please refer to the `Numpy documentation <http://www.scipy.org/Tentative_NumPy_Tutorial#head-864862d3f2bb4c32f04260fac61eb4ef34788c4c>`.
You can also check on the shape of the array with ::

    my_map.shape

If you'd like to rip the data out of the SunPy map to use elsewhere
you can use::

    import numpy as np
    var = np.array(aia_map)
    
The variable var is now simply a `Numpy ndarray <http://docs.scipy.org/doc/numpy/reference/arrays.ndarray.html>`. Fortunately this is not generally necessary
as the SunPy map inherits directly from the NumPy ndarray and therefore applying NumPy methods
directly on the SunPy map is generally safe. So one could, for example, apply the following code::

    import numpy as np
    log_map = np.log(my_map)

The np.log() function will be applied on each element in the map data array as expected. A few NumPy
ndarray methods are also directly available such as, for example, to find the maximum and mean 
in the map data just use::

    my_map.min()
    my_map.mean()

This is the advantage of inheritance!

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
The make_map() method described above can also handle a list of maps. If the maps are
of a different type (e.g. different instruments) than the result of make_map is 
what we call a Composite Map. So for example to create a simple composite map::

    my_maps = sunpy.make_map([sunpy.EIT_195_IMAGE, sunpy.RHESSI_IMAGE])

A Composite map is different from a regular map and therefore different associated methods.
To list which maps are part of your composite map use::

    my_maps.list_maps()

Similar to all SunPy data objects, the composite map also has an associated show() method and a 
number of associated methods to customize your plot. For example, the following code turns 
adds a new map, sets its transparency to 25%, turns on contours from 50% to 90% for the second map, 
and then plots the result::

    my_maps.add_map(sunpy.AIA_171_IMAGE)
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

