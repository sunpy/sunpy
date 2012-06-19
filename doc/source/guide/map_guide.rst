----
Maps
----

Maps in SunPy are spatially-aware data arrays. In other words, they are 2-dimensional data associated with a coordinate system. In this guide, we will 
cover some of the basic functionality of maps. Once you've read through this guide check out the :doc:`code reference</reference/index>` for a more
thorough look at SunPy maps.

1. Creating maps
----------------
SunPy contains a number of example fits files. To make things easy,
SunPy includes several example files which are used throughout the docs. These
files have names like `sunpy.AIA_171_IMAGE` and `sunpy.RHESSI_IMAGE`.
To create the sample AIA map type the following into your interactive Python shell::

	import sunpy
	aia_map = sunpy.make_map(sunpy.AIA_171_IMAGE)

The variable aia is a SunPy map. To create a map from a local fits try
something like the following ::

    my_map = sunpy.make_map('/mydirectory/mymap.fits')

SunPy should automatically detect the type of file (e.g. fits), what instrument it is 
associated with (e.g. AIA, EIT, LASCO) and will automatically look in the appropriate places for the fits
keywords it needs to interpret the coordinate system. If the type of fits file 
is not recognized then SunPy will try some default fits keywords but results
may vary. SunPy can also create maps from the jpg2000 files from
`helioviewer.org <http://helioviewer.org/>`.

2. Inspecting maps
------------------
A map contains a number of data-associated attributes. To get of your map summary simply
type::

    my_map
    
This will show a representation of the data as well as its associated
attributes. A number of attributes are also available, for example the date, 
exposure time, map center, xrange, yrange
other::

    my_map.date
    my_map.exptime
    my_map.center
    my_map.xrange
    my_map.yrange
    
To get a list of all of the attributes check the documentation by typing::

	help(my_map)
	
The header using the function::

    my_map.get_header()
    
This returns a dictionary with the header information as read from the source
file. You can also easily check the yrange and xrange with the following 
commands::

3. Getting at the data
----------------------
As a SunPy map inherits from the NumPy ndarray it can be accessed in the same
way through indexing. To get the 0th element in the array ::

    aia_map[0,0]
    aia_map[0][0]
    
For more information about this please refer to the `Numpy documentation < http://www.scipy.org/Tentative_NumPy_Tutorial#head-864862d3f2bb4c32f04260fac61eb4ef34788c4c>`.

If you'd like to rip the data (ndarray) out of the SunPy map to use elsewhere
you can use::

    import numpy as np
    var = np.array(aia_map)
    
The variable var is now a NumPy ndarray. 
The SunPy map inherits from the NumPy ndarray and therefore also has access to all of
the ndarray functionality allowing further inspection and transformation of the 
data. So for example to find the maximum and minimum in the map data just use::

    my_map.min()
    my_map.max()

4. Creating a plot of your map
------------------------------
The SunPy map object has its own built-in plot function so that it is easy to
quickly view your map. To create a plot just type::

	my_map.show()
	
In addition, to enable users to modify the map it is also possible to grab the
matplotlib figure object by using the plot() command instead of the show() 
command. This makes it possible to use the SunPy plot as the foundation for a 
more complicated figure. For example,


5. Working with your map
------------------------
Part of the philosophy of the map object is to provide most of the basic
functionality that a scientist would want therefore a map also contains a number
of map-specific methods such as resizing a map or grabbing a subview. To get 
a list of the methods available for a map type::

	help(my_map)
	
and check out the methods section!

