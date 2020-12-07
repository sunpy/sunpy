"""
===============================================
Creating an Offset frame using NorthOffsetFrame
===============================================

This is an example to show some possible usage of ``NorthOffsetFrame``.
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
# Creating the offset frame  
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
