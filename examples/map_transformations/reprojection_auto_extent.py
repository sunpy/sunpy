"""
==================================
Reprojecting with Automatic Extent
==================================

This example demonstrates how you can automatically determine the extent
needed for the output of a reprojection instead of specifying it ahead of time.
"""
# sphinx_gallery_thumbnail_number = 3

import matplotlib.pyplot as plt

import astropy.units as u
from astropy.coordinates import SkyCoord

import sunpy.map
from sunpy.coordinates import SphericalScreen
from sunpy.data.sample import AIA_171_IMAGE

######################################################################
# Let's build off of the final reprojection made in the example
# :ref:`sphx_glr_generated_gallery_map_transformations_reprojection_spherical_screen.py`.

aia_map = sunpy.map.Map(AIA_171_IMAGE)
aia_map.plot_settings['norm'].vmin = 0
aia_map.plot_settings['norm'].vmax = 10000

new_observer = SkyCoord(70*u.deg, 20*u.deg, 1*u.AU, obstime=aia_map.date,
                        frame='heliographic_stonyhurst')

out_shape = aia_map.data.shape

out_ref_coord = SkyCoord(0*u.arcsec, 0*u.arcsec, obstime=new_observer.obstime,
                         frame='helioprojective', observer=new_observer,
                         rsun=aia_map.coordinate_frame.rsun)
out_header = sunpy.map.make_fitswcs_header(
    out_shape,
    out_ref_coord,
    scale=u.Quantity(aia_map.scale),
    instrument=aia_map.instrument,
    wavelength=aia_map.wavelength
)

with SphericalScreen(aia_map.observer_coordinate, only_off_disk=True):
    outmap_screen_off_disk = aia_map.reproject_to(out_header)

fig = plt.figure()
ax = fig.add_subplot(projection=outmap_screen_off_disk)
outmap_screen_off_disk.plot(axes=ax)

######################################################################
# Note that some of the data is clipped at the top and bottom because
# ``out_shape`` was not defined large enough in the vertical
# direction. On the other hand, there is too much empty space to the
# left and right because ``out_shape`` was defined too large in the
# horizontal direction, which actually results in unnecessary
# computation time. That is also why one should avoid specifying
# a conservatively large ``out_shape``.
#
# Instead, we can enable the automatic determination of extent by
# :meth:`~sunpy.map.GenericMap.reproject_to` via the keyword
# ``auto_extent``.  Let's first try ``auto_extent='edges'``, which
# is often sufficient for many situations. This setting will
# precompute the reprojection of just the edges of the map to see
# what the output extent should be.

with SphericalScreen(aia_map.observer_coordinate, only_off_disk=True):
    outmap_edges = aia_map.reproject_to(out_header, auto_extent='edges')

fig = plt.figure()
ax = fig.add_subplot(projection=outmap_edges)
outmap_edges.plot(axes=ax)

######################################################################
# In this case, ``auto_extent='edges'`` misses the fact that the Sun
# in the middle part of the map actually extends to the left of one of
# edges in the reprojection. Thus, instead, we need to specify
# ``auto_extent='all'`` to precompute the reprojection of all pixels
# in the map.  This setting takes more time, but will ensure that no
# part of the reprojected map is clipped.

with SphericalScreen(aia_map.observer_coordinate, only_off_disk=True):
    outmap_all = aia_map.reproject_to(out_header, auto_extent='all')

fig = plt.figure()
ax = fig.add_subplot(projection=outmap_all)
outmap_all.plot(axes=ax)

plt.show()
