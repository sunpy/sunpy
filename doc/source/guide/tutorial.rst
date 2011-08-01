----------------
A Quick Tutorial
----------------

A simple example of how to plot a map::

	import sunpy
	import matplotlib.cm as cm
	import matplotlib.colors as colors
	map = sunpy.Map(“doc/sample-data/AIA20110319_105400_0171.fits”)
	map.plot(cmap=cm.hot, norm=colors.Normalize(1, 2048))

and there you go!