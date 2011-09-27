----------------
A Quick Tutorial
----------------

A simple example of how to plot a map::

	import sunpy
	from matplotlib import cm
	from matplotlib import colors
	aia = sunpy.Map(sunpy.AIA_171_IMAGE)
	aia.show(cmap=cm.hot, norm=colors.Normalize(1, 2048))

and there you go!

Colormaps
---------

There are a number of color maps defined in SunPy which are used for data from particular missions (e.g. SDO/AIA). 
A simple example on how to use the color maps provided by SunPy: ::

	import sunpy.cm as cm
	# cmlist is a dictionary with all of the color tables
	# to list all of the keys of the dictionary
	cm.cmlist.keys()
	# to grab a particular colortable then
	cmap = cm.cmlist.get('sdoaia94')
	# you can also get a visual representation of all of the color tables 
	cm.show_colormaps()

Solar Physical Constants
------------------------

SunPy contains a convienient list of solar-related physical constants. Here is a short bit of code to
get you started: ::
	
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

Anytim and Utilities
--------------------
SunPy also contains a number of utility functions which may be useful in general. To grab those import them as such: ::
	from sunpy.util import util as util
For example, parsing a standard time strings: ::
	util.anytim('2004/02/05 12:00')
This returns a datetime object. All SunPy functions which require time as an input sanitize the input using util.anytim. To get the day of year: ::
	util.day_of_year('2004-Jul-05 12:00:02')
or the julian day: ::
	util.julian_day((2010,4,30))
