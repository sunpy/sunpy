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
