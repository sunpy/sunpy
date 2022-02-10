"""
=========================================
Overplotting HMI Contours on an AIA Image
=========================================

This example shows how to use `~astropy.visualization.wcsaxes` to overplot
unaligned HMI magnetic field strength contours on an AIA map.
"""
# sphinx_gallery_thumbnail_number = 2

import matplotlib.pyplot as plt
import numpy as np

import astropy.units as u
from astropy.coordinates import SkyCoord

import sunpy.map
from sunpy.data.sample import AIA_193_IMAGE, HMI_LOS_IMAGE

################################################################################
# First let's load two of the sample files into two Map objects.

aia, hmi = sunpy.map.Map(AIA_193_IMAGE, HMI_LOS_IMAGE)

################################################################################
# To make the plot neater, we start by submapping the same region.
# We define the region in HGS coordinates and then apply the same submap to
# both the HMI and AIA maps.

bottom_left = SkyCoord(30 * u.deg, -40 * u.deg, frame='heliographic_stonyhurst')
top_right = SkyCoord(70 * u.deg, 0 * u.deg, frame='heliographic_stonyhurst')

sub_aia = aia.submap(bottom_left, top_right=top_right)
sub_hmi = hmi.submap(bottom_left, top_right=top_right)

################################################################################
# To highlight the fact that the AIA and HMI images are not aligned, let us
# quickly view the two maps side-by-side.

fig = plt.figure(figsize=(11, 5))

ax1 = fig.add_subplot(121, projection=sub_aia)
sub_aia.plot(axes=ax1, clip_interval=(1, 99.99)*u.percent)

ax2 = fig.add_subplot(122, projection=sub_hmi)
sub_hmi.plot(axes=ax2)

################################################################################
# In the next plot we will start by plotting the same aia submap, and draw a
# heliographic grid on top.

# sphinx_gallery_defer_figures

fig = plt.figure(figsize=(8, 8))

ax = fig.add_subplot(projection=sub_aia)
sub_aia.plot(axes=ax, clip_interval=(1, 99.99)*u.percent)
sub_aia.draw_grid(axes=ax)

ax.set_title("AIA 193 with HMI magnetic field strength contours", y=1.1)

################################################################################
# Now we want to draw the contours, to enhance the appearance of the plot we
# explicitly list the levels, but then make them symmetric around 0

# sphinx_gallery_defer_figures

levels = [50, 100, 150, 300, 500, 1000] * u.Gauss

################################################################################
# matplotlib requires the levels to be sorted, so we order them from lowest to
# highest by reversing the array.

# sphinx_gallery_defer_figures

levels = np.concatenate((-1 * levels[::-1], levels))

################################################################################
# Before we add the contours to the axis we store the existing bounds of the
# as overplotting the contours will sometimes change the bounds, we re-apply
# them to the axis after the contours have been added.

# sphinx_gallery_defer_figures

bounds = ax.axis()

################################################################################
# We use the map method `~.GenericMap.draw_contours` to simplify this process,
# but this is a wrapper around `~matplotlib.pyplot.contour`. We set the
# colormap, line width and transparency of the lines to improve the final
# appearance.

# sphinx_gallery_defer_figures

cset = sub_hmi.draw_contours(levels, axes=ax, cmap='seismic', alpha=0.5)
ax.axis(bounds)

################################################################################
# Finally, add a colorbar. We add an extra tick to the colorbar at the 0 point
# to make it clearer that it is symmetric, and tweak the size and location to
# fit with the axis better.

plt.colorbar(cset,
             label=f"Magnetic Field Strength [{sub_hmi.unit}]",
             ticks=list(levels.value) + [0],
             shrink=0.8,
             pad=0.17)

plt.show()
