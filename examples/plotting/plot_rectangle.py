"""
==========================================================
Drawing a Rectangle on a map
==========================================================

In this example we are going to look at how we can plot rectangles on maps
while describing the dimensions using different methods.
"""

################################################################################
# Fundamentally, we will be using the `~sunpy.map.GenericMap.draw_rectangle` method to
# draw the rectangles. Its function signature is like so:
#
# `draw_rectangle(bottom_left, *, width: Unit(‘deg’) = None,
# height: Unit(‘deg’) = None, axes=None, top_right=None, **kwargs)`
#
# All of the extra kwargs are passed on to the `~matplotlib.patches.Rectangle` instance,
# which allows us to modify the color, linestyle, etc. of the plotted rectangle.
# A full overview of available options can be found in the Rectangle documentation.
#
# Therefore, the different forms in which an input can be provided to the function are:
#
# * `~astropy.coordinates.SkyCoord` ((x1, x2), (y1, y2))
# * `~astropy.coordinates.SkyCoord` bottom_left, `~astropy.coordinates.SkyCoord` top_right
# * `~astropy.coordinates.SkyCoord` bottom_left, `~astropy.units.Quantity` width, `~astropy.units.Quantity` height
#
# Note that from SunPy v2.1, every argument to `~sunpy.map.GenericMap.draw_rectangle` apart from
# 'bottom_left' (i.e. 'top_right', 'width' and 'height') should be passed as a keyword argument.
#
# We also have a fourth option - to specify the dimensions in terms of pixel coordinates.
#
# To do this, we need to first create a `~sunpy.map.GenericMap.submap` by passing
# `~astropy.units.Quantity` objects with the units as pixel. Then we can use the
# `~sunpy.map.GenericMap.bottom_left_coord` and `~sunpy.map.GenericMap.top_right_coord` properties
# as parameters for `~sunpy.map.GenericMap.draw_rectangle`.

import matplotlib.pyplot as plt

import astropy.units as u
from astropy.coordinates import SkyCoord

import sunpy.data.sample
import sunpy.map

aia_1 = sunpy.map.Map(sunpy.data.sample.AIA_171_IMAGE)
aia_2 = sunpy.map.Map(sunpy.data.sample.AIA_171_IMAGE)
aia_3 = sunpy.map.Map(sunpy.data.sample.AIA_171_IMAGE)
aia_4 = sunpy.map.Map(sunpy.data.sample.AIA_171_IMAGE)

fig = plt.figure(figsize=(12, 12))

# ((x1, x2), (y1, y2))
ax_1 = fig.add_subplot(221, projection=aia_1)
coords = SkyCoord((-400, 200)*u.arcsec, (-300, 100)*u.arcsec, frame=aia_1.coordinate_frame)
aia_1.draw_rectangle(coords, axes=ax_1, color='blue', linestyle='-')
aia_1.plot()

# bottom_left, top_right
ax_2 = fig.add_subplot(222, projection=aia_2)
bottom_left = SkyCoord(-400*u.arcsec, -300*u.arcsec, frame=aia_2.coordinate_frame)
top_right = SkyCoord(300*u.arcsec, 200*u.arcsec, frame=aia_2.coordinate_frame)
aia_2.draw_rectangle(bottom_left, top_right=top_right, axes=ax_2, color='green', linestyle='--')
aia_1.plot()

# bottom_left, width, height
ax_3 = fig.add_subplot(223, projection=aia_3)
bottom_left = SkyCoord(-500*u.arcsec, -200*u.arcsec, frame=aia_3.coordinate_frame)
width = 700 * u.arcsec
height = 500 * u.arcsec
aia_3.draw_rectangle(bottom_left, width=width, height=height, axes=ax_3, color='yellow', linestyle='-.')
aia_3.plot()

# pixel coordinates
ax_4 = fig.add_subplot(224, projection=aia_4)
bottom_left_pixel = [400, 400] * u.pixel
top_right_pixel = [700, 600] * u.pixel
aia_submap = aia_4.submap(bottom_left_pixel, top_right=top_right_pixel)
aia_4.draw_rectangle(aia_submap.bottom_left_coord, top_right=aia_submap.top_right_coord, axes=ax_4, color='red', linestyle=':')
aia_4.plot()

plt.show()
