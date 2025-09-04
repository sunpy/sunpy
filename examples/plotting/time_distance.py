"""
=====================================================
Creating a time-distance plot from a sequence of maps
=====================================================

This example showcases how you can use :func:`sunpy.map.pixelate_coord_path`
and :func:`sunpy.map.sample_at_coords` on a sequence of images to create a
time-distance diagram accounting for solar differential rotation using
`sunpy.coordinates.propagate_with_solar_surface` and dealing with off-disk
pixels using `sunpy.coordinates.screens.SphericalScreen`
"""
# sphinx_gallery_thumbnail_number = -1
import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import numpy as np

import astropy.units as u
from astropy.coordinates import SkyCoord

from sunpy.coordinates import propagate_with_solar_surface
from sunpy.coordinates.screens import SphericalScreen
from sunpy.map import Map, pixelate_coord_path, sample_at_coords
from sunpy.net import Fido
from sunpy.net import attrs as a

###############################################################################
# First, we will need to acquire a sequence of images from SDO/AIA.
# We will use a data from 2012 containing a nice example of loop oscillations.

# To keep the online build and download low using a short time range expand to
# 18:05 - 18:30 to see more of the oscillation.

query = Fido.search(
    a.Time('2012-10-20 18:14:00', '2012-10-20 18:19:00'),
    a.Instrument.aia,
    a.Wavelength(171*u.angstrom),
)
files = Fido.fetch(query, site="NSO")
files = sorted(files)


###############################################################################
# Our target will be a set of loops in the corona above an active region.
# First load the FITS files we downloaded above in to map sequence and create a
# rectangular submap or cutout around the area of interest. We will also define
# the path, in this case a line, along which we want to make the time-distance
# plot.

aia_seq = Map(files)
corner = SkyCoord(Tx=-1150*u.arcsec, Ty=-500*u.arcsec,
                  frame=aia_seq[0].coordinate_frame)
ref_sub_map = aia_seq[0].submap(corner, width=250*u.arcsec, height=450*u.arcsec)

line_coords = SkyCoord([-1030, -1057]*u.arcsec, [-220, -206]*u.arcsec,
                       frame=aia_seq[0].coordinate_frame)

###############################################################################
# Next we can plot the full disk map, over-plotting the region of interest and
# also plot the submap and the line along which we want to make the
# time-distance plot.

fig = plt.figure(figsize=(8, 5))
full_disk_ax = fig.add_subplot(121, projection=aia_seq[0])
aia_seq[0].plot(axes=full_disk_ax, clip_interval=[1,99]*u.percent)
aia_seq[0].draw_quadrangle(corner, width=250*u.arcsec, height=450*u.arcsec,
                           axes=full_disk_ax)
sub_map_ax = fig.add_subplot(122, projection=ref_sub_map)
ref_sub_map.plot(axes=sub_map_ax, clip_interval=[1,99]*u.percent)
sub_map_ax.plot_coord(line_coords)

###############################################################################
# There are two approaches that can be used to extract time distance
# measurements from the data:
#
# 1. Reproject all the maps to a common world coordinate system (WCS)
# 2. Transform the coordinates of the line extracted from the reference map
#    to the coordinate systems of the subsequent maps
#
# We will use both and show they achieve almost the same results choosing which
# approach to use will depend on the exact use case.
#
# We will start of with the first approach and reproject all the maps to common WCS
# of the ``ref_sub_map.wcs`` while also taking account of differential rotation of
# using `~sunpy.coordinates.propagate_with_solar_surface` and off-disk pixels using
# `~sunpy.coordinates.screens.SphericalScreen`.

reprojected_sub_maps = []
for cur_map in aia_seq:
    cur_map = cur_map/cur_map.exposure_time
    with (propagate_with_solar_surface(),
          SphericalScreen(cur_map.observer_coordinate, only_off_disk=True)):
        reprojected_sub_maps.append(cur_map.reproject_to(ref_sub_map.wcs))
reprojected_sub_maps = Map(reprojected_sub_maps, sequence=True)

###############################################################################
# Now that we have reprojected all the maps to common WCS we can extract the
# pixel coordinates once using :func:`~sunpy.map.pixelate_coord_path` to
# determine the coordinates for every pixel that intersects with the physical
# path and then use :func:`~sunpy.map.sample_at_coords` sample the data at
# these coordinates.

# As the maps are all aligned only need to extract the coordinates once
intensity_coords = pixelate_coord_path(aia_seq[0], line_coords)

intensities_reproject = []

for aia_map in reprojected_sub_maps:
    intensities_reproject.append(sample_at_coords(aia_map, intensity_coords).value)

###############################################################################
# For the second approach we need to transform the ``intensity_coords`` into
# the coordinate frame of each map and then extract the data at corresponding
# pixel coordinates.

intensities_transform = []
for cur_map in aia_seq:
    with propagate_with_solar_surface(), SphericalScreen(cur_map.observer_coordinate,
                                                         only_off_disk=True):
        # The coordinate will automatically be transformed into the cur_map frame
        intensities_transform.append(
            (sample_at_coords(cur_map, intensity_coords)/cur_map.exposure_time).value)

###############################################################################
# Now we have obtained the raw data we need to prepare it for platting and
# calculate the extents of the x and y acies for the final plot.

# The above will give us a list of 1D arrays, one for each map in the sequence.
# We need to stack them into a 2D array.
intensities_reproject = np.vstack(intensities_reproject)
intensities_transform = np.vstack(intensities_transform)

# This defines the distance along the path in arcseconds.
angular_separation = intensity_coords.separation(intensity_coords[0]).to(u.arcsec)

# Get correct values for the extent
extent = [aia_seq[0].date.datetime, aia_seq[-1].date.datetime,
          0, angular_separation[-1].value]

###############################################################################
# Plot the reference submap, line and extracted data from both methods and
# the difference between them.

fig = plt.figure(figsize=(10, 5), layout="constrained")
left, right = fig.subfigures(nrows=1, ncols=2, width_ratios=[0.6, 0.75])
left_ax = left.add_subplot(111, projection=reprojected_sub_maps[0])
right_ax = right.subplot_mosaic([['repro'], ['trans'], ['diff']],
                                sharex=True, sharey=True)

imag_ax = reprojected_sub_maps[0].plot(axes=left_ax, clip_interval=[1, 99]*u.percent)
plt.colorbar(imag_ax, ax=left_ax)
left_ax.plot_coord(line_coords)

right_ax['repro'].imshow(
    intensities_reproject.T, aspect='auto', interpolation='none',
    extent=extent, cmap=imag_ax.get_cmap()
)
plt.colorbar(right_ax['repro'].images[0], ax=right_ax['repro'])
right_ax['trans'].imshow(
    intensities_transform.T, aspect='auto', interpolation='none',
    extent=extent, cmap=imag_ax.get_cmap()
)
plt.colorbar(right_ax['trans'].images[0], ax=right_ax['trans'])
right_ax['diff'].imshow(
    (intensities_reproject-intensities_transform).T, interpolation='none',
    aspect='auto', extent=extent, cmap='bwr'
)
plt.colorbar(right_ax['diff'].images[0], ax=right_ax['diff'])

locator = mdates.AutoDateLocator(minticks=4)
formatter = mdates.ConciseDateFormatter(locator)
right_ax['diff'].xaxis.set_major_locator(locator)
right_ax['diff'].xaxis.set_major_formatter(formatter)
right_ax['diff'].xaxis.set_minor_locator(mdates.MinuteLocator())

right_ax['repro'].set_title('Reproject')
right_ax['trans'].set_title('Transform')
right_ax['diff'].set_title('Difference')

right_ax['diff'].set_xlabel('Time [UTC]')
right_ax['trans'].set_ylabel('Distance [arcsec]')

plt.show()
