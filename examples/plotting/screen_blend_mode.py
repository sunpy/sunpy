"""
======================
Blending maps together
======================

This example shows how to blend two maps.

``matplotlib`` by itself provides only alpha-based transparency for
superimposing one image onto another, which may not be visually appealing.
For better results, one can manipulate the rendered image data arrays to perform a
standard image-compositing `blend mode <https://en.wikipedia.org/wiki/Blend_modes>`__.
"""
import matplotlib.pyplot as plt

import astropy.units as u

import sunpy.data.sample
import sunpy.map
from sunpy.coordinates import SphericalScreen

###############################################################################
# Let's load two maps for blending. We reproject the second map to the
# coordinate frame of the first map for proper compositing, taking care to use
# the :class:`~sunpy.coordinates.SphericalScreen` context manager in order to
# preserve off-disk data.

a171 = sunpy.map.Map(sunpy.data.sample.AIA_171_IMAGE)
a131 = sunpy.map.Map(sunpy.data.sample.AIA_131_IMAGE)
with SphericalScreen(a171.observer_coordinate):
    a131 = a131.reproject_to(a171.wcs)

###############################################################################
# Let's first plot the two maps individually to check the desired plot
# settings.

fig1 = plt.figure(figsize=(10, 4))
ax1 = fig1.add_subplot(121, projection=a171)
ax2 = fig1.add_subplot(122, projection=a131)

a171.plot(axes=ax1, clip_interval=(1, 99.999)*u.percent)
a131.plot(axes=ax2, clip_interval=(1, 99.9)*u.percent)

###############################################################################
# Now let's prepare for the computation of the blended image. We plot both
# maps using the desired plot settings in order to obtain the ``matplotlib``
# image instances for future use.

fig2 = plt.figure()
ax = fig2.add_subplot(projection=a171)

im171 = a171.plot(axes=ax, clip_interval=(1, 99.99)*u.percent)
im131 = a131.plot(axes=ax, clip_interval=(1, 99.9)*u.percent)

# sphinx_gallery_defer_figures

###############################################################################
# We use those ``matplotlib`` image instances to obtain RGB data arrays after
# normalization (including clipping) and the color maps have been applied. We
# need to provide the current ``matplotlib`` renderer and turn off resampling
# to screen pixels. To facilitate later arithmetic, we then scale the data
# arrays from being integers ranging 0-255 to being floats ranging 0.0-1.0.

renderer = fig2.canvas.get_renderer()
rgb171 = im171.make_image(renderer, unsampled=True)[0] / 255.
rgb131 = im131.make_image(renderer, unsampled=True)[0] / 255.

# sphinx_gallery_defer_figures

###############################################################################
# With these RGB arrays, we can now perform the calculation for the
# `screen blend mode <https://en.wikipedia.org/wiki/Blend_modes#Screen>`__.
# A different blend mode would simply require a different expression.

rgb_composite = 1 - (1 - rgb171) * (1 - rgb131)

# sphinx_gallery_defer_figures

###############################################################################
# Finally, we plot the composite image on the same axes.

ax.imshow(rgb_composite, origin='lower')
ax.set_title('Composite using the screen blend mode')
plt.show()

# sphinx_gallery_thumbnail_number = 2
