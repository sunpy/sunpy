"""
==============================================
Creating a visualization with ImageAnimatorWCS
==============================================

This example shows how to create a simple visualization using
`~sunpy.visualization.animator.ImageAnimatorWCS`.
"""
# Start by importing the necessary modules.
from itertools import product

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

# This dictionary comprehension is extracting the three basic keywords we need
# to create a astropy.wcs.WCS header: 'CTYPE','CUNIT' and 'CDELT'
# from the meta information stored in the 'map_sequence'.
wcs_input_dict = {f'{key}{n+1}': map_sequence.all_meta()[0].get(f'{key}{n}')
                  for n, key in product([1, 2], ['CTYPE', 'CUNIT', 'CDELT'])}

# Now we need to get the time difference between the two observations.
t0, t1 = map(parse_time, [k['date-obs'] for k in map_sequence.all_meta()])
time_diff = (t1 - t0).to(u.s)
wcs_input_dict.update({'CTYPE1': 'TIME', 'CUNIT1': time_diff.unit.name, 'CDELT1': time_diff.value})

# We can now just pass this into astropy.wcs.WCS to create our WCS header.
wcs = astropy.wcs.WCS(wcs_input_dict)

# Now the resulting WCS object will look like:
print(wcs)

###############################################################################
# Now we can create the animation.
# `~sunpy.visualization.animator.ImageAnimatorWCS` assumes the last two axes
# are the two that form the image. However, they are the first two in this case,
# so we change this by passing in the ``image_axes`` keyword.

wcs_anim = ImageAnimatorWCS(sequence_array, wcs=wcs, vmax=1000, image_axes=[0, 1])

###############################################################################
# You might notice that the animation could do with having the axes look neater.
# As we are using a WCS object, we can use `~astropy.visualization.wcsaxes`
# to change the properties of the visualization. The documentation is
# `here <https://docs.astropy.org/en/stable/visualization/wcsaxes/index.html>`_
# and details each possible method. In this example, we will just showcase some
# of these possible methods.

# We have to recreate the visualization since we displayed it earlier.
wcs_anim = ImageAnimatorWCS(sequence_array, wcs=wcs, vmax=1000, image_axes=[0, 1])

# We need to access the underlying WCSAxes object, so we can change the
# properties of our visualization. The methods on the ``axes.coords`` attribute
# of ``wcs_anim`` allows us to set various properties.
# We assign them to a easy to remember variable name.
# Note the order of assignment out from ``wcs_anim.axes.coords``.

time, solar_y, solar_x = wcs_anim.axes.coords

# Now we can label the X and Y axes.
solar_x.set_axislabel('Solar X (arsec)')
solar_y.set_axislabel('Solar Y (arsec)')

# Move the axis labels to avoid the slider.
solar_x.set_axislabel_position('t')
solar_y.set_axislabel_position('r')

# Now we can change the spacing and we have to use `astropy.units` here.
# We are setting the spacing to be quite small here and this will cause overlap.
solar_x.set_ticks(spacing=1*u.arcmin, color='black')
solar_y.set_ticks(spacing=1*u.arcmin, color='black')

# We can make sure the ticks do not overlap.
solar_x.set_ticklabel(exclude_overlapping=True)
# Set this one to false, so we can see the difference.
solar_y.set_ticklabel(exclude_overlapping=False)
