----------------
A Quick Tutorial
----------------

A simple example of how to plot a map::

	import sunpy
	import matplotlib.cm as cm
	import matplotlib.colors as colors
	map = sunpy.Map(sunpy.AIA_171_IMAGE)
	map.show(cmap=cm.hot, norm=colors.Normalize(1, 2048))

and there you go!

How to make an image from a RHESSI calibrated event list::

	from pylab import *
	import sunpy
	from sunpy.instr import rhessi
	result = rhessi.backprojection(sunpy.RHESSI_EVENT_LIST)
	imshow(result)
