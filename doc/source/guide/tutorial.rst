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
