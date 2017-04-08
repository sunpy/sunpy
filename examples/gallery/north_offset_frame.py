"""
===============
North Offset Frame
===============

It is used to create a frame derived from Heliographic, which has the north pole at some point of interest. 
In this frame, lines of longitude form circles radially away from the point, and lines of latitude measure angular distance from the point.
"""

import sunpy.coordinates as sc
import astropy.units as u
from astropy.coordinates import SkyCoord

import sunpy.coordinates
from sunpy.coordinates.offset_frame import NorthOffsetFrame
import sunpy.map
from sunpy.data.sample import AIA_171_ROLL_IMAGE, AIA_171_IMAGE
from astropy.wcs import WCS
import astropy.units as u
%matplotlib notebook
import matplotlib.pyplot as plt
import numpy as np

north = SkyCoord(20*u.deg, 20*u.deg, frame="heliographic_stonyhurst")
f = NorthOffsetFrame(north=north)
m = sunpy.map.Map(AIA_171_ROLL_IMAGE)

north = SkyCoord(20*u.deg, 20*u.deg, frame="heliographic_stonyhurst")
f = NorthOffsetFrame(north=north)

fig = plt.figure(figsize=(10,10))
ax = fig.add_subplot(111, projection=m)
m.plot()

overlay = ax.get_coords_overlay(f)
overlay[0].set_ticks(spacing=20. * u.deg, color='white')
overlay.grid(ls='-', color='red')
###################################################
# In this example the new frame is shifted so the new north pole is at (20, 20) in the Heliographic Stonyhurst frame. The new grid is overplotted in red.



#####################################################################################
# Now let's try and select a range of pixels based on longitude in this frame, we want a raidal wedge of pixels away from our new north pole.
# Get arrays of all pixel coordinates for the whole image.

x, y = np.meshgrid(*[np.arange(v.value) for v in m.dimensions])*u.pix

#####################################################################################
# Calculate HPC coords for all pixels

rot = SkyCoord(*m.pixel_to_data(x,y), frame=m.coordinate_frame)

#####################################################################################
# Transform to our shifted frame

rot = rot.transform_to(f)

####################################################################################
# Now select a band of pixels based on

seg = np.logical_and(rot.lon<10*u.deg, rot.lon>0*u.deg).nonzero()
m.data[seg] = 0
m.peek()
