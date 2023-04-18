"""
=======================================================
Segmenting a Map based on transformation of coordinates
=======================================================

This example demonstrates extracting a region of a particular map based on
world coordinates in different systems.
"""
# sphinx_gallery_thumbnail_number = 2
import matplotlib.pyplot as plt
import numpy as np

import astropy.units as u
from astropy.coordinates import SkyCoord

import sunpy.map
from sunpy.coordinates import NorthOffsetFrame
from sunpy.data.sample import AIA_171_IMAGE

######################################################################
# We start with the sample data.

smap = sunpy.map.Map(AIA_171_IMAGE)

######################################################################
# A utility function gives us access to the helioprojective coordinate of each
# pixel in this map. From this we then transform the coordinates to the required frame.
# For this example we are going to extract a region based on the
# heliographic Stonyhurst coordinates, so we transform to that frame.

all_hpc = sunpy.map.all_coordinates_from_map(smap)
all_hgs = all_hpc.transform_to("heliographic_stonyhurst")

######################################################################
# Let's then segment the data based on coordinates and create a boolean mask
# where `True` indicates invalid or deselected data.
# Numpy's masked arrays allow for a combination of standard numpy array and
# a mask array. When an element of the mask is `True`, the corresponding element
# of the associated array is said to be masked (invalid).
# For more information about numpy's masked arrays see :mod:`numpy.ma`.
# We now mask out all values not in our coordinate range or where the
# coordinates are NaN (because they could not be transformed to the
# surface of the Sun).

segment_mask = np.logical_or(all_hgs.lon >= 35 * u.deg, all_hgs.lon <= -35 * u.deg)
segment_mask |= np.isnan(all_hgs.lon)

######################################################################
# To plot the segment separately, we create a new map with the segment as the mask.

new_frame_map = sunpy.map.Map(smap.data, smap.meta, mask=segment_mask)
fig = plt.figure()
ax = fig.add_subplot(projection=new_frame_map)
new_frame_map.plot(axes=ax)
new_frame_map.draw_grid(axes=ax, color='red')
plt.show()

######################################################################
# We can perform various mathematical operations on the extracted segment such
# as averaging the pixel values or finding the sum of the segment.

masked_data = np.ma.array(new_frame_map.data, mask=new_frame_map.mask)

print(f"Original Map : mean = {smap.data.mean()}, sum = {smap.data.sum()}")
print(f"Segment : mean = {masked_data.mean()}, sum = {masked_data.sum()}")

######################################################################
# Using `sunpy.coordinates.NorthOffsetFrame`
# ------------------------------------------
# Let us offset the north pole and create the frame.

north = SkyCoord(20 * u.deg, 20 * u.deg, frame="heliographic_stonyhurst")
offset_frame = NorthOffsetFrame(north=north)

######################################################################
# We then transform coordinates to the offsetted frame and segment the data
# based on conditions.

all_hpc = sunpy.map.all_coordinates_from_map(smap)
offsetted_coords = all_hpc.transform_to(offset_frame)
segment_mask = np.logical_or(offsetted_coords.lon >= 30 * u.deg,
                             offsetted_coords.lon <= -20 * u.deg)

######################################################################
# Masking out the NaN values of ``offsetted_coords.lon``, we get:

segment_mask |= np.isnan(offsetted_coords.lon)

######################################################################
# Let's plot the offsetted segment separately.

offsetted_map = sunpy.map.Map(smap.data, smap.meta, mask=segment_mask)
fig = plt.figure()
ax = fig.add_subplot(projection=smap)
offsetted_map.plot(axes=ax)
overlay = ax.get_coords_overlay(offset_frame)
overlay[0].set_ticks(spacing=30. * u.deg)
overlay.grid(ls='--', color='blue')
offsetted_map.draw_grid(axes=ax, color='red')
plt.show()

######################################################################
# We can also find the maximum, minimum or average pixel values of the segment
# and compare it with the original map.

offset_masked_data = np.ma.array(offsetted_map.data, mask=offsetted_map.mask)
print(f"Original Map : mean = {smap.data.mean()}, "
      f"maximum value = {smap.data.max()}, "
      f"minimum value = {smap.data.min()}")
print(f"Offset segment : mean = {offset_masked_data.mean()}, "
      f"maximum value = {offset_masked_data.max()}, "
      f"minimum value = {offset_masked_data.min()}")
