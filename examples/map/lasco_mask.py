"""
=================================
Creating a mask for LASCO C2 data
=================================

In this example, we will manually create a mask to block the occulter in an
unprocessed LASCO C2 coronagraph image.
"""
import matplotlib.pyplot as plt
import numpy as np

import astropy.units as u
from astropy.coordinates import SkyCoord

from sunpy.map import Map
from sunpy.map.maputils import all_coordinates_from_map
from sunpy.net import Fido
from sunpy.net import attrs as a

###############################################################################
# First, download some unprocessed LASCO C2 data with `~sunpy.net.Fido`.

result = Fido.search(a.Time('2011/06/07 06:30', '2011/06/07 06:36'),
                     a.Instrument.lasco,
                     a.Detector.c2)
lasco_file = Fido.fetch(result)
lasco_map = Map(lasco_file)

###############################################################################
# The LASCO C2 coronagraph has a field of view extending from 2-6 solar
# radii. So, our mask will have two parts: an inner component which masks the
# occulter and an outer component which masks data outside the field of view.
#
# We will follow a process similar to
# :ref:`sphx_glr_generated_gallery_computer_vision_techniques_finding_masking_bright_pixels.py`
# to express the coordinates relative to the occulter center.

pixel_coords = all_coordinates_from_map(lasco_map)
solar_center = SkyCoord(0*u.deg, 0*u.deg, frame=lasco_map.coordinate_frame)
pixel_radii = np.sqrt((pixel_coords.Tx-solar_center.Tx)**2 +
                      (pixel_coords.Ty-solar_center.Ty)**2)
# Note that the inner mask extends just beyond 2 solar radii to mask the
# Fresnel diffraction caused by the occulter edge.
mask_inner = pixel_radii < lasco_map.rsun_obs*2.4
mask_outer = pixel_radii > lasco_map.rsun_obs*6
final_mask = mask_inner + mask_outer

###############################################################################
# To apply the final mask, we must create a new map.

masked_lasco = Map(lasco_map.data, lasco_map.meta, mask=final_mask)

# Before plotting the map, we need to create a new colormap to ensure we mask
# the bad values correctly.
occult_colormap = lasco_map.cmap.copy()
occult_colormap.set_bad('black')

fig = plt.figure()
ax1 = fig.add_subplot(1, 2, 1, projection=lasco_map)
ax2 = fig.add_subplot(1, 2, 2, projection=masked_lasco)

lasco_map.plot(clip_interval=(2, 98)*u.percent, cmap=occult_colormap, axes=ax1)
lasco_map.draw_limb()
ax1.set_title("Level 1 LASCO C2")

masked_lasco.plot(clip_interval=(2, 98)*u.percent, cmap=occult_colormap, axes=ax2)
masked_lasco.draw_limb()
ax2.set_title("Masked LASCO C2")

plt.show()
