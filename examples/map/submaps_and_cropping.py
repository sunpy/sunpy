"""
=================
Creating a submap
=================

How to cut a rectangular region from a map with
:meth:`~sunpy.map.GenericMap.submap`.

:meth:`~sunpy.map.GenericMap.submap` is not the only way to cut down a map.
:meth:`~ndcube.NDCube.crop` does a related job and is shown in
:ref:`sphx_glr_generated_gallery_map_cropping_a_map.py`; for which to reach for,
see :ref:`sunpy-topic-guide-submap-vs-crop`.
"""
# sphinx_gallery_tags = ["Map", "SWAP", "Coordinates"]

import matplotlib.pyplot as plt

import astropy.units as u
from astropy.coordinates import SkyCoord

import sunpy.data.sample
import sunpy.map

################################################################################
# We start with the sample data

swap_map = sunpy.map.Map(sunpy.data.sample.SWAP_LEVEL1_IMAGE)

###############################################################################
# To cut out a rectangular region of a Map, create a submap by specifying its
# top-right and bottom-left corners as `~astropy.coordinates.SkyCoord` objects.

top_right = SkyCoord(0 * u.arcsec, -200 * u.arcsec, frame=swap_map.coordinate_frame)
bottom_left = SkyCoord(-900 * u.arcsec, -900 * u.arcsec, frame=swap_map.coordinate_frame)
swap_submap = swap_map.submap(bottom_left, top_right=top_right)

###############################################################################
# Let's plot the results.

fig = plt.figure()
ax = fig.add_subplot(projection=swap_submap)
image = swap_submap.plot(axes=ax)
swap_submap.draw_limb(axes=ax)
swap_submap.draw_grid(axes=ax)

# Make some room and put the title at the top of the figure
ax.set_position([0.1, 0.1, 0.8, 0.7])
ax.set_title(ax.get_title(), pad=45)

plt.show()
