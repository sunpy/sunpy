----------------
A Quick Tutorial
----------------

Welcome to the SunPy tutorial! This brief tutorial will walk you through some 
of the functionality currently offered by SunPy. Start by reading this tutorial
and trying out some of the examples demonstrated. Once you've completed the
tutorial check out the :doc:`code reference</reference/index>` for a more
thorough look at the functionality available.

1. Plotting
-----------

Let's begin by creating a simple plot of an AIA image. To make things easy,
SunPy includes several example files which are used throughout the docs. These
files have names like `sunpy.AIA_171_IMAGE` and `sunpy.RHESSI_IMAGE`.

Try typing the below example into your interactive Python shell::

	import sunpy
	from matplotlib import cm
	from matplotlib import colors
	aia = sunpy.Map(sunpy.AIA_171_IMAGE)
	aia.show(cmap=cm.hot, norm=colors.Normalize(1, 2048))

If everything has been configured properly you should see a standard-looking
AIA 171 image with a colorbar on the right-hand side and a title and some 
labels.

There is lot going on here, but we will walk you through the example. Briefly,
in the first few lines we are just importing SunPy and a couple other plotting
related modules that we will use in the example. On the fourth line we create a
SunPy Map object which is basically just a spatially-aware image or data array.
On the last line we then plot the map object, adding a couple additional
parameters to specify a color map to use and how we wish to scale the image.

Over the next few sections we will explain some of these features in more depth
and then move onto more other modules included in SunPy.

Specifying a Colormap
^^^^^^^^^^^^^^^^^^^^^

There are a number of color maps defined in SunPy which are used for data from 
particular missions (e.g. SDO/AIA). 
A simple example on how to use the color maps provided by SunPy: ::

	import sunpy.cm as cm
	
	# cmlist is a dictionary with all of the color tables
	# to list all of the keys of the dictionary
	cm.cmlist.keys()

	# to grab a particular colortable then
	cmap = cm.cmlist.get('sdoaia94')

	# you can also get a visual representation of all of the color tables 
	cm.show_colormaps()

2. Solar Physical Constants
---------------------------

SunPy contains a convienient list of solar-related physical constants. Here is 
a short bit of code to get you started: ::
	
	from sunpy.sun import constants as con

	# one astronomical unit (the average distance between the Sun and Earth)
	print(con.au)

	# the solar radius
	print(con.radius)

	# not all constants have a shortcut assigned to them (as above)
	# the rest of the constants are stored in a dictionary
	solar_constants = con.physical_constants

	# to get a list of all of the values stored in this dictionary
	solar_constants.keys()
	
	# or you can use the following convinience method to list them all
	con.print_all()

3. Anytim and Utilities
-----------------------

SunPy also contains a number of utility functions which may be useful in 
general. Here is a short example: ::

	from sunpy.util import util as util
	
	# parsing a standard time strings
	util.anytim('2004/02/05 12:00')
	
	# This returns a datetime object. All SunPy functions which require 
	# time as an input sanitize the input using util.anytim. 	
	util.day_of_year('2004-Jul-05 12:00:02')
	
	# the julian day
	util.julian_day((2010,4,30))
	
4. Querying the VSO
-------------------
There are a couple different ways to query and download data from the VSO using
SunPy. The method you should use depends first on your preference with respect
to query style: the main method of querying uses a syntax that is unique to
SunPy and may require some getting used to, but is extremely flexible and
powerful. To make it easy for people coming from SSW to get started, a second
"legacy" API also exists which works is very much the same way as VSO_GET in
IDL.

Further, for each of the two query APIs there are interactive and
non-interactive versions available, depending on the type of work you are doing.

The below example demonstrates a simple query for STEREO EUVI data using the
non-interactive version of the main API::

    from sunpy.net import vso
    
    # create a new VSOClient instance
    client = vso.VSOClient()
    
    # build our query
    result = client.query(
        vso.attrs.Time((2011, 9, 20), (2011, 9, 21)),
        vso.attrs.Instrument('euvi')
    )
    
    # print the number of matches
    print("Number of records found: %d " % result.no_records())
   
    # download matches to /download/path
    res = client.get(result, path="/download/path").wait()

Note that specifying a path is optional and if you do not specify one the files
will simply be downloaded into a temporary directory (e.g. /tmp/xyz).

