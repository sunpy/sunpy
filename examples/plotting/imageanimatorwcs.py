"""
==============================================
Creating a visualization with ImageAnimatorWCS
==============================================

This example shows how to create a simple visualization using
`~sunpy.visualization.animator.ImageAnimatorWCS`.
"""
# Start by importing the necessary modules.
from itertools import product

import matplotlib.pyplot as plt

import astropy.wcs
import astropy.units as u

import sunpy.map
from sunpy.time import parse_time
from sunpy.visualization.animator import ImageAnimatorWCS
from sunpy.data.sample import AIA_171_IMAGE, AIA_193_IMAGE

################################################################################
# To showcase how to visualize a sequence of 2D images using
# `~sunpy.visualization.animator.ImageAnimatorWCS`, we will use images from
# our sample data. The problem with this is that they are not part of
# a continuous dataset. To overcome this we wil do two things.
# Create a stacked array of the images and create a `~astropy.wcs.WCS` header.
# The easiest method  for the array is to create a `~sunpy.map.MapSequence`.

# Here we only use two files but you could pass in a larger selection of files.
map_sequence = sunpy.map.Map(AIA_171_IMAGE, AIA_193_IMAGE, sequence=True)
# Now we can just cast the sequence away into a NumPy array.
sequence_array = map_sequence.as_array()

###############################################################################
# Now we need to create the `~astropy.wcs.WCS` header that
# `~sunpy.visualization.animator.ImageAnimatorWCS` will need.
# To create the new header we can use the stored meta information from the
# ``map_sequence``.
# We will be using a overally complex method to create a header.

# This dictionary comphersion is extracting the three basic keywords we need
# to create a astropy.wcs.WCS header: 'CTYPE','CUNIT', 'CDELT', 'CRVAL'
# and 'CRPIX' from the meta information stored in the 'map_sequence'.
wcs_input_dict = {f'{key}{n+1}': map_sequence.all_meta()[0].get(f'{key}{n}')
                  for n, key in product([1, 2], ['CTYPE', 'CUNIT', 'CDELT', 'CRVAL', 'CRPIX'])}

# Now we need to get the time difference between the two observations.
t0, t1 = map(parse_time, [k['date-obs'] for k in map_sequence.all_meta()])
time_diff = (t1 - t0).to(u.s)
wcs_input_dict.update({'CTYPE1': 'Time', 'CUNIT1': time_diff.unit.name, 'CDELT1': time_diff.value})

# Now add the correct number of axes to the WCS header
wcs_input_dict.update({'NAXIS': 3, 'NAXIS1': 2, 'NAXIS2': 1024, 'NAXIS3': 1024})

# We can now just pass this into astropy.wcs.WCS to create our WCS header.
wcs = astropy.wcs.WCS(wcs_input_dict)

# Now the resulting WCS object will look like:
print(wcs)

###############################################################################
# Now we can create the animation.
# As `~sunpy.visualization.animator.ImageAnimatorWCS` assumes the last two axes
# are the two form the image, we can change this  by passing in the
# ``image_axes`` keyword.

#sequence_array = np.rollaxis(sequence_array, -1)
wcs_anim = ImageAnimatorWCS(sequence_array, wcs=wcs, vmax=1000, image_axes=[0, 1])

###############################################################################
# You might notice that the animation could do with having the axes ticks
# being tidied up. To do this we can call
# `~sunpy.visualization.animator.ImageAnimatorWCS.update_axes_style`.

wcs_anim = ImageAnimatorWCS(sequence_array, wcs=wcs, vmax=1000, image_axes=[0, 1])
wcs_anim.update_axes_style((20 * u.arcsec), 't', 'r')
# Notice the title for this one.
# Is this the bug for https://github.com/sunpy/sunpy/pull/2894?

###############################################################################
# You might notice that the animation could do with having the axes ticks
# being tidied up. Since we are using a WCS object, we can use
# `~astropy.visualization.wcsaxes` to change the properties of the visualization.
# `Documentation is here <http://docs.astropy.org/en/stable/visualization/wcsaxes/index.html>_`

wcs_anim = ImageAnimatorWCS(sequence_array, wcs=wcs, vmax=1000, image_axes=[0, 1])

# We need to access the underlying WCS axes object.
# This way we change various properties of our visualization.
ax = wcs_anim.axes

# Now we can label the X and Y axes.
ax.set_xlabel('Solar X (arsec)')
ax.set_ylabel('Solar Y (arsec)')

# To change the spacing of the ticks for each axis we have to go deeper.
# ax.coords holds each axis in a list and we can index.
y_axis, x_axis = ax.coords[1], ax.coords[2]

# Move the axis labels to avoid the slider.
y_axis.set_axislabel_position('r')
x_axis.set_axislabel_position('t')

# Now we can change the spacing.
y_axis.set_ticks(spacing=1*u.arcmin, color='black')
x_axis.set_ticks(spacing=1*u.arcmin, color='black')

# We can make sure the ticks do not overlap.
y_axis.set_ticklabel(exclude_overlapping=True)
# Set this one to false, so we can see the difference.
x_axis.set_ticklabel(exclude_overlapping=True)
