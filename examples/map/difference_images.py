"""
===========================
Plotting a difference image
===========================

This example shows how to compute and plot a difference image.
"""
# sphinx_gallery_thumbnail_number = 2
import matplotlib.colors as colors
import matplotlib.pyplot as plt

from astropy.visualization import ImageNormalize, SqrtStretch

import sunpy.data.sample
import sunpy.map

###########################################################################
# When analyzing solar imaging data it is often useful to look at the
# difference from one time step to the next (running difference) or the
# difference from the start of the sequence (base difference). In this
# example, we'll use a sequence of AIA 193 cutout images taken during a
# flare.
#
# First, load a series of images into a `~sunpy.map.MapSequence`.

m_seq = sunpy.map.Map([
    sunpy.data.sample.AIA_193_CUTOUT01_IMAGE,
    sunpy.data.sample.AIA_193_CUTOUT02_IMAGE,
    sunpy.data.sample.AIA_193_CUTOUT03_IMAGE,
    sunpy.data.sample.AIA_193_CUTOUT04_IMAGE,
    sunpy.data.sample.AIA_193_CUTOUT05_IMAGE,
], sequence=True)

###########################################################################
# Let's take a look at each image in the sequence. Note that these images
# are sampled at a 6.4 minute cadence, much lower than the actual 12 s
# AIA cadence. We first adjust the plot settings on each image to ensure
# the colorbar is the same at each time step.
for m in m_seq:
    m.plot_settings['norm'] = ImageNormalize(vmin=0, vmax=5e3, stretch=SqrtStretch())

plt.figure()
ani = m_seq.plot()
plt.show()

###########################################################################
# And now we can take the actual difference. We will compute the running
# difference and the base difference for each image in the sequence. For the
# case of the first entry in the running difference, we just subtract it
# from itself.
#
# But we have to decide what to do with the metadata. For example, what
# time does this difference image correspond to? The time of the first or
# second image? The mean time? You'll have to decide what makes most sense
# for your application. Here, by subtracting the data of the first (and
# previous) data array from each map in the sequence, the resulting
# difference maps have the same metadata as each corresponding map in the
# sequence.
#
# Note that, because arithmetic operations between `~sunpy.map.GenericMap`
# objects are not supported, we subtract just the array data (with units
# attached) of the second map from the first map. The ``quantity`` attribute
# returns the image data as an `~astropy.units.Quantity`, where the resulting
# units are those returned by the ``unit`` attribute of the map.

m_seq_base = sunpy.map.Map([m - m_seq[0].quantity for m in m_seq[1:]], sequence=True)
m_seq_running = sunpy.map.Map(
    [m - prev_m.quantity for m, prev_m in zip(m_seq[1:], m_seq[:-1])],
    sequence=True
)

###########################################################################
# Finally, let's plot the difference maps.
# We'll apply a colormap and re-normalize the intensity so that it shows
# up well.
# First, we show the base difference map.

norm = colors.Normalize(vmin=-200, vmax=200)
for m in m_seq_base:
    m.plot_settings['norm'] = norm
    m.plot_settings['cmap'] = 'Greys_r'
plt.figure()
ani = m_seq_base.plot(title='Base Difference')
plt.colorbar(extend='both', label=m_seq_base[0].unit.to_string())
plt.show()

###########################################################################
# Then, we show the running difference map.

for m in m_seq_running:
    m.plot_settings['norm'] = norm
    m.plot_settings['cmap'] = 'Greys_r'
plt.figure()
ani = m_seq_running.plot(title='Running Difference')
plt.colorbar(extend='both', label=m_seq_running[0].unit.to_string())
plt.show()
