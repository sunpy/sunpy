"""
=======================================================
Segmenting a Map based on transformation of coordinates
=======================================================

A demostration creating a segment of a particular map from
transformed coordinates.
"""
import matplotlib.pyplot as plt
import numpy as np
from reproject import reproject_interp

import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.wcs import WCS

import sunpy.map
from sunpy.coordinates import Helioprojective, NorthOffsetFrame, RotatedSunFrame, transform_with_sun_center
from sunpy.data.sample import AIA_171_IMAGE

###############################################################################
# We start with the sample data.
smap = sunpy.map.Map(AIA_171_IMAGE)

###############################################################################
# A utility function gives us access to the helioprojective coordinate of each
# pixel. From this we then transform the coordinates to the required frame.
all_hpc = sunpy.map.all_coordinates_from_map(smap)
all_hgs = all_hpc.transform_to("heliographic_stonyhurst")

###############################################################################
# Let's then segment the data based on coordinates and create a boolean mask
# where True indicates the segmented mask.
segment = np.logical_and(all_hgs.lon < 10 * u.deg, all_hgs.lon > 0 * u.deg)

###############################################################################
# Let's plot the segment separately, we create a new map with the segment as the
# mask.
new_frame_map = sunpy.map.Map(
    smap.data,
    smap.meta,
    mask=np.logical_not(segment))
fig = plt.figure()
fig.add_subplot(projection=smap)
new_frame_map.plot()
plt.show()

###############################################################################
# We can perform various mathematical operations on the extracted segment such
# as averaging the pixel values or finding the sum of the segment.
mean_pix = np.average(segment)
sum_pix = np.sum(segment)

###############################################################################
# This example can also be extended while using NorthOffsetFrame
# Let us offset the north pole and create the frame
north = SkyCoord(20 * u.deg, 20 * u.deg, frame="heliographic_stonyhurst")
f = NorthOffsetFrame(north=north)

###############################################################################
# We then transform coordinates to the Offsetted frame and segment the data
# based on conditions
all_hpc = sunpy.map.all_coordinates_from_map(smap)
ofsetted_coords = all_hpc.transform_to(f)
segment = np.logical_and(
    ofsetted_coords.lon < 30 * u.deg,
    ofsetted_coords.lon > -20 * u.deg)

###############################################################################
# Let's plot the ofsetted segment separately
ofsetted_map = sunpy.map.Map(
    smap.data,
    smap.meta,
    mask=np.logical_not(segment))
fig = plt.figure()
fig.add_subplot(projection=smap)
ofsetted_map.plot()
plt.show()

###############################################################################
# We can also find the maximum, minimum or average pixel values of the segment
mean_pix = np.average(segment)
max_pix = np.max(segment*smap.data)
min_pix = np.min(segment*smap.data)
