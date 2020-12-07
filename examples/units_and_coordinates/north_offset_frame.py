"""
===============================================
Creating an Offset frame using NorthOffsetFrame
===============================================

This example shows a possible use of `~sunpy.coordinates.NorthOffsetFrame`

A common use for this is to create a frame derived from Heliographic, which
has the north pole at some point of interest. In this new frame, lines of
longitude form great circles radially away from the point, and lines of
latitude measure angular distance from the point.
In this example the new frame is shifted so the new north pole is at (20,
20) in the Heliographic Stonyhurst frame. The new grid is overplotted in
blue.

"""

##############################################################################
# Import the required modules.
import matplotlib.pyplot as plt
import astropy.units as u

from astropy.coordinates import SkyCoord
from sunpy.coordinates import NorthOffsetFrame

import sunpy.map
from sunpy.data.sample import AIA_171_IMAGE

##############################################################################
# Now to create the offset frame. 
m = sunpy.map.Map(AIA_171_IMAGE)

north = SkyCoord(20*u.deg, 20*u.deg, frame="heliographic_stonyhurst")
new_frame = NorthOffsetFrame(north=north)

ax = plt.subplot(projection=m)
m.plot()

overlay = ax.get_coords_overlay('heliographic_stonyhurst')
overlay[0].set_ticks(spacing=30. * u.deg, color='white')
overlay.grid(ls='-', color='white')

overlay = ax.get_coords_overlay(new_frame)
overlay[0].set_ticks(spacing=30. * u.deg)
overlay.grid(ls='-', color='blue')
