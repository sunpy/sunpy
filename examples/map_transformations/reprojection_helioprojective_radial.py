"""
============================================================
Reprojecting a Helioprojective Map to Helioprojective Radial
============================================================

In this example we use `reproject` to transform a helioprojective
`~sunpy.map.Map` to helioprojective radial.
"""
# sphinx_gallery_thumbnail_number = 2

import matplotlib.pyplot as plt

import astropy.units as u

import sunpy.data.sample
import sunpy.map
from sunpy.coordinates import HelioprojectiveRadial
from sunpy.map.header_helper import make_hpr_header
from sunpy.visualization import show_hpr_impact_angle

###############################################################################
# Let's start with an AIA map in the helioprojective Cartesian
# (`~sunpy.coordinates.Helioprojective`) frame.

aia_map = sunpy.map.Map(sunpy.data.sample.AIA_193_IMAGE)

###############################################################################
# We now plot the map and overlay grid lines for position angle in the
# helioprojective radial (`~sunpy.coordinates.HelioprojectiveRadial`) frame.

fig = plt.figure()
ax = fig.add_subplot(projection=aia_map)
aia_map.plot(axes=ax)

overlay = ax.get_coords_overlay(
    HelioprojectiveRadial(
        obstime=aia_map.reference_date,
        observer=aia_map.observer_coordinate
    )
)

overlay['delta'].add_tickable_gridline('limb', 960*u.arcsec - 90*u.deg)

overlay['psi'].grid(color='red', linestyle='dashed')
overlay['psi'].set_ticks(spacing=30*u.deg)
overlay['psi'].set_ticks_visible(False)
overlay['psi'].set_ticklabel(color='red')
overlay['psi'].set_ticklabel_position(('limb',))
overlay['psi'].set_axislabel('Helioprojective Position Angle', color='red')
overlay['psi'].set_axislabel_position('r')

overlay['delta'].set_ticks_visible(False)
overlay['delta'].set_ticklabel_visible(False)
overlay['delta'].set_axislabel('')

###############################################################################
# In order to reproject this map to a helioprojective radial frame, we use
# :func:`sunpy.map.make_hpr_header` to create an appropriate FITS WCS header.

hpr_header = make_hpr_header(
    aia_map.observer_coordinate,
    (1000, 360),
    1.8*u.arcsec
)

###############################################################################
# With the new header, reproject the data into the new coordinate frame.
# The :meth:`~sunpy.map.GenericMap.reproject_to` defaults to using
# the fast :func:`reproject.reproject_interp` algorithm, but a different
# algorithm can be specified (e.g., :func:`reproject.reproject_adaptive`).

outmap = aia_map.reproject_to(hpr_header)

###############################################################################
# Now we plot the reprojected map. First we add an overlay for the original
# helioprojective Cartesian frame, but disable its associated ticks and labels.

fig = plt.figure()
ax = fig.add_subplot(projection=outmap)
outmap.plot(axes=ax)

ax.coords[0].set_ticks(spacing=30*u.deg)
ax.set_aspect(0.3)

overlay = ax.get_coords_overlay(aia_map.coordinate_frame)
overlay.grid(color='blue', linestyle='dashed')
for coord in overlay:
    coord.set_ticks(spacing=500*u.arcsec)
    coord.set_ticks_visible(False)
    coord.set_ticklabel_visible(False)
    coord.set_axislabel('')

# sphinx_gallery_defer_figures

###############################################################################
# Due to the FITS WCS machinery, the impact-angle component of helioprojective
# radial is instead stored as a declination by subtracting 90 degrees, and the
# vertical axis being declination can be unfamiliar. To change the declination
# tick labels to impact-angle tick labels, we use
# :func:`~sunpy.visualization.show_hpr_impact_angle` to modify
# the rendering of the tick labels on the appropriate axis (``ax.coords[1]``).

show_hpr_impact_angle(ax.coords[1])
ax.coords[1].set_format_unit('arcsec')

plt.show()
