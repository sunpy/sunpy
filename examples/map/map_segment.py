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
from sunpy.data.sample import AIA_171_ROLL_IMAGE

###############################################################################
# We start with the sample data.
smap = sunpy.map.Map(AIA_171_ROLL_IMAGE)

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
# Another example would be to transform the coordinates to a RotatedSunFrame
# Referring to the example of differentially rotating a map, we create a
# RotatedSunFrame to transform the coordinates
in_time = smap.date

out_time = in_time + 5 * u.day
out_frame = Helioprojective(observer='earth', obstime=out_time)
rot_frame = RotatedSunFrame(base=out_frame, rotated_time=in_time)
out_shape = smap.data.shape
out_center = SkyCoord(0 * u.arcsec, 0 * u.arcsec, frame=out_frame)
header = sunpy.map.make_fitswcs_header(out_shape,
                                       out_center,
                                       scale=u.Quantity(smap.scale))
out_wcs = WCS(header)
out_wcs.coordinate_frame = rot_frame

with transform_with_sun_center():
    arr, _ = reproject_interp(smap, out_wcs, out_shape)

out_warp = sunpy.map.Map(arr, out_wcs)
out_warp.plot_settings = smap.plot_settings

###############################################################################
# Transforming the coordinates to the RotatedSunFrame
all_hpc = sunpy.map.all_coordinates_from_map(out_warp)
rot_coords = all_hpc.transform_to(rot_frame)

###############################################################################
# Creating the segment and converting it into a separate map
segment = np.logical_and(
    rot_coords.Tx > 800 * u.arcsec,
    rot_coords.Ty < 300 * u.arcsec)

rotated_map = sunpy.map.Map(
    out_warp.data,
    out_warp.meta,
    mask=np.logical_not(segment))
rotated_map.plot_settings = out_warp.plot_settings

###############################################################################
# Let's plot the new segment
fig = plt.figure()
fig.add_subplot(projection=smap)
rotated_map.plot()
plt.show()
