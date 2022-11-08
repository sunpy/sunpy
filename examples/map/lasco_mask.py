"""
=================================
Creating a mask for LASCO C2 data
=================================

In this example, we will manually create a mask to block the occulter in an
unprocessed LASCO C2 coronagraph image.
"""
import matplotlib.pyplot as plt
import numpy as np
import numpy.ma as ma

import astropy.units as u

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
# First, let's define the Sun center using ``CRPIXi`` and apply a shift (if
# necessary) to the reference coordinates ``CRVALi``
# (see `here <https://lasco-www.nrl.navy.mil/level_1/level_1_keywords.html>`_).

x_occult = lasco_map.meta['crval1']/lasco_map.meta['cdelt1'] + lasco_map.meta['crpix1']
y_occult = lasco_map.meta['crval2']/lasco_map.meta['cdelt2'] + lasco_map.meta['crpix2']
occult_coord = lasco_map.pixel_to_world(x_occult*u.pixel, y_occult*u.pixel)

###############################################################################
# Next, we will follow a process similar to
# :ref:`sphx_glr_generated_gallery_computer_vision_techniques_finding_masking_bright_pixels.py`
# to express the coordinates relative to the occulter center. Then, we will
# use `numpy.ma` functions to mask the correct regions.
# Note that the inner mask extends just beyond 2 solar radii to mask the
# Fresnel diffraction caused by the occulter edge.

hpc_coords = all_coordinates_from_map(lasco_map)
norm_mask = np.sqrt((hpc_coords.Tx-occult_coord.Tx)**2 +
                    (hpc_coords.Ty-occult_coord.Ty)**2)
mask_inner = ma.masked_less_equal(norm_mask, lasco_map.rsun_obs*2.4)
mask_outer = ma.masked_greater_equal(norm_mask, lasco_map.rsun_obs*6)
final_mask = mask_inner + mask_outer

###############################################################################
# To apply the final mask, create a new map using ``.mask`` to get the actual
# masked array.

masked_lasco = Map(lasco_map.data, lasco_map.meta, mask=final_mask.mask)

# Before plotting the map, we need to create a new colormap to ensure we mask
# the bad values correctly.
occult_colormap = lasco_map.cmap.copy()
occult_colormap.set_bad('black')

plt.figure()
masked_lasco.plot(clip_interval=(2, 98)*u.percent, cmap=occult_colormap)
masked_lasco.draw_limb()
plt.title("Custom mask for LASCO C2")

plt.show()
