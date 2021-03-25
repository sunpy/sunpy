"""
=======================================================
Segmenting a Map based on transformation of coordinates
=======================================================

A demonstration extracting a region of a particular map based on
world coordinates in different systems.
"""
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
# pixel. From this we then transform the coordinates to the required frame.
# For this example we are going to extract a region based on the
# heliographic Stonyhurst coordinates, so we transform to that frame.
all_hpc = sunpy.map.all_coordinates_from_map(smap)
all_hgs = all_hpc.transform_to("heliographic_stonyhurst")

######################################################################
# Let's then segment the data based on coordinates and create a boolean mask
# where True indicates invalid or deselected data.
# This is the convention used by the `numpy.ma` module.
segment = np.logical_or(
    all_hgs.lon >= 50 * u.deg,
    all_hgs.lon <= 0 * u.deg)

######################################################################
# Let's plot the segment separately, we create a new map with the segment as the
# mask.
# Numpy's masked arrays allow for a combination of standart numpy array and
# a mask array. When an element of the mask is True, the corresponding element
# of the associated array is said to be masked (invalid).
# More information about Numpy's masked arrays :
# https://numpy.org/doc/stable/reference/maskedarray.generic.html#using-numpy-ma
# From the segment, we now mask out all values where `all_hgs.lon` is
# NaN.
segment = np.logical_or(
    segment, np.isnan(
        all_hgs.lon))

new_frame_map = sunpy.map.Map(
    smap.data,
    smap.meta,
    mask=segment)
fig = plt.figure()
ax = plt.subplot(projection=smap)
new_frame_map.plot()
new_frame_map.draw_grid()
plt.show()

######################################################################
# We can perform various mathematical operations on the extracted segment such
# as averaging the pixel values or finding the sum of the segment.
new_frame_mask = np.ma.array(smap.data, mask=segment)

print(
    f"Original Map : mean = {smap.data.mean()}, sum = {smap.data.sum()}")
print(
    f"Segmented Map : mean = {new_frame_mask.mean()}, sum = {new_frame_mask.sum()}")

######################################################################
# This example can also be extended while using NorthOffsetFrame
# Let us offset the north pole and create the frame
north = SkyCoord(
    20 * u.deg,
    20 * u.deg,
    frame="heliographic_stonyhurst")
f = NorthOffsetFrame(north=north)

######################################################################
# We then transform coordinates to the Offsetted frame and segment the data
# based on conditions
all_hpc = sunpy.map.all_coordinates_from_map(smap)
ofsetted_coords = all_hpc.transform_to(f)
segment = np.logical_or(
    ofsetted_coords.lon >= 30 * u.deg,
    ofsetted_coords.lon <= -20 * u.deg)

######################################################################
# Masking out the NaN values of `ofsetted_coords.lon`, we get
segment = np.logical_or(
    segment, np.isnan(
        ofsetted_coords.lon))

######################################################################
# Let's plot the ofsetted segment separately
offsetted_map = sunpy.map.Map(
    smap.data,
    smap.meta,
    mask=segment)
fig = plt.figure()
ax = plt.subplot(projection=smap)
offsetted_map.plot()
overlay = ax.get_coords_overlay(f)
overlay[0].set_ticks(spacing=30. * u.deg)
overlay.grid(ls='--', color='blue')
offsetted_map.draw_grid()
plt.show()

######################################################################
# We can also find the maximum, minimum or average pixel values of the segment
# and compare it with the original map.
offset_mask = np.ma.array(smap.data, mask=segment)
print(f"Original Map : mean = {smap.data.mean()}, "
      f"maximum value = {smap.data.max()}, "
      f"minimum value = {smap.data.min()}")
print(f"Segmented Map : mean = {offset_mask.mean()}, "
      f"maximum value = {offset_mask.max()}, "
      f"minimum value = {offset_mask.min()}")
