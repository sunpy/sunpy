"""
==============================================
Creating a visualization with ArrayAnimatorWCS
==============================================

This example shows how to create a simple visualization using
`~sunpy.visualization.animator.ArrayAnimatorWCS`.
"""
# Start by importing the necessary modules.

import matplotlib.pyplot as plt

import astropy.units as u
import astropy.wcs

import sunpy.map
from sunpy.data.sample import AIA_171_IMAGE, AIA_193_IMAGE
from sunpy.time import parse_time
from sunpy.visualization.animator import ArrayAnimatorWCS

################################################################################
# To showcase how to visualize a sequence of 2D images using
# `~sunpy.visualization.animator.ArrayAnimatorWCS`, we will use images from
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
# `~sunpy.visualization.animator.ArrayAnimatorWCS` will need.
# To create the new header we can use the stored meta information from the
# ``map_sequence``.

# Now we need to get the time difference between the two observations.
t0, t1 = map(parse_time, [k['date-obs'] for k in map_sequence.all_meta()])
time_diff = (t1 - t0).to(u.s)

m = map_sequence[0]

wcs = astropy.wcs.WCS(naxis=3)
wcs.wcs.crpix = u.Quantity([0*u.pix] + list(m.reference_pixel))
wcs.wcs.cdelt = [time_diff.value] + list(u.Quantity(m.scale).value)
wcs.wcs.crval = [0, m._reference_longitude.value, m._reference_latitude.value]
wcs.wcs.ctype = ['TIME'] + list(m.coordinate_system)
wcs.wcs.cunit = ['s'] + list(m.spatial_units)
wcs.wcs.aux.rsun_ref = m.rsun_meters.to_value(u.m)

# Now the resulting WCS object will look like:
print(wcs)

###############################################################################
# Now we can create the animation.
# `~sunpy.visualization.animator.ArrayAnimatorWCS` requires you to select which
# axes you want to plot on the image. All other axes should have a ``0`` and
# sliders will be created to control the value for this axis.

wcs_anim = ArrayAnimatorWCS(sequence_array, wcs, [0, 'x', 'y'], vmax=1000)

plt.show()

###############################################################################
# You might notice that the animation could do with having the axes look
# neater. `~sunpy.visualization.ArrayAnimatorWCS` provides a way of setting
# some display properties of the `~astropy.visualization.wcsaxes.WCSAxes`
# object on every frame of the animation via use of the ``coord_params`` dict.
# They keys of the ``coord_params`` dict are either the first half of the
# ``CTYPE`` key, the whole ``CTYPE`` key or the entries in
# ``wcs.world_axis_physical_types`` here we use the short ctype identifiers for
# the latitude and longitude axes.


coord_params = {
    'hpln': {
        'axislabel': 'Helioprojective Longitude',
        'ticks': {'spacing': 10*u.arcmin, 'color': 'black'}
    },
    'hplt': {
        'axislabel': 'Helioprojective Latitude',
        'ticks': {'spacing': 10*u.arcmin, 'color': 'black'}
    },
}


# We have to recreate the visualization since we displayed it earlier.
wcs_anim = ArrayAnimatorWCS(sequence_array, wcs, [0, 'x', 'y'], vmax=1000,
                            coord_params=coord_params)

plt.show()
