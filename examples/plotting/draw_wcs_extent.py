"""
===========================
Drawing the Extent of a WCS
===========================

This example demonstrates how to draw the extent of a WCS on a `~sunpy.map.Map`
"""
# sphinx_gallery_thumbnail_number = 3
import matplotlib.pyplot as plt

import astropy.units as u
from astropy.coordinates import SkyCoord

import sunpy.data.sample
import sunpy.map
import sunpy.visualization.drawing

################################################################################
# First, load an AIA map and create a submap around a band of active regions.

m_aia = sunpy.map.Map(sunpy.data.sample.AIA_193_JUN2012)
m_aia_sub = m_aia.submap(SkyCoord(Tx=-500*u.arcsec, Ty=-500*u.arcsec, frame=m_aia.coordinate_frame),
                         top_right=SkyCoord(Tx=700*u.arcsec, Ty=0*u.arcsec, frame=m_aia.coordinate_frame))

fig = plt.figure()
ax = fig.add_subplot(projection=m_aia_sub)
m_aia_sub.plot(axes=ax)

################################################################################
# To show the context of the submap, we can draw the extent of the submap on top
# of the full-disk map.

fig = plt.figure()
ax = fig.add_subplot(projection=m_aia)
m_aia.plot(axes=ax)
m_aia_sub.draw_extent(axes=ax)

################################################################################
# Additionally, we can draw the extent of the submap on top of an EUV
# observation from STEREO A which was separated from SDO by over 116 degrees.
# Note that using `sunpy.visualization.drawing.extent` allows for drawing the
# extent of any WCS which contains two celestial axes.

m_euvi = sunpy.map.Map(sunpy.data.sample.STEREO_A_195_JUN2012)
fig = plt.figure()
ax = fig.add_subplot(projection=m_euvi)
m_euvi.plot(axes=ax)
sunpy.visualization.drawing.extent(ax, m_aia_sub.wcs)

plt.show()
