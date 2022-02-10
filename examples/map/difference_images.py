"""
===========================
Plotting a difference image
===========================

This example shows how to compute and plot a difference image.
"""
# sphinx_gallery_thumbnail_number = 2
import matplotlib.colors as colors
import matplotlib.pyplot as plt

import astropy.units as u

import sunpy.data.sample
import sunpy.map

###########################################################################
# When analyzing solar imaging data, particularly for flares, it is often
# useful to look at the difference from one time step to the next (running
# difference) or the difference from the start of the sequence (base
# difference). In this example, we'll use a sequence of AIA 193 cutout
# images taken during a flare.
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
# AIA cadence.

fig = plt.figure(figsize=(12, 3))
for i, m in enumerate(m_seq):
    ax = fig.add_subplot(1, 5, i+1, projection=m)
    delta_t = (m.date - m_seq[0].date).to('minute')
    m.plot(axes=ax, clip_interval=(0.5, 99.9)*u.percent,
           title=f'{delta_t:.2f}')
    ax.coords[0].set_axislabel('Solar-X' if i == 0 else ' ')
    ax.coords[0].set_ticklabel_visible(i == 0)
    ax.coords[1].set_axislabel('Solar-Y' if i == 0 else ' ')
    ax.coords[1].set_ticklabel_visible(i == 0)
plt.show()

###########################################################################
# And now we can take the actual difference. We will compute the running
# difference and the base difference for the last image in the sequence.
#
# But we have to decide what to do with the metadata. For example, what
# time does this difference image correspond to? The time of the first or
# second image? The mean time? You'll have to decide what makes most sense
# for your application. Here, by subtracting the data of the first (and
# second to last) data array from the last map in the sequence, the resulting
# difference map has the same metadata as the last map in the sequence.

m_diff_base = m_seq[-1] - m_seq[0].quantity
m_diff_running = m_seq[-1] - m_seq[-2].quantity

###########################################################################
# Finally, let's plot the difference maps.
# We'll apply a colormap and re-normalize the intensity so that it shows
# up well.
# First, we show the base difference map.

norm = colors.Normalize(vmin=-200, vmax=200)
plt.figure()
m_diff_base.plot(cmap='Greys_r', norm=norm, title='Base Difference')
plt.colorbar(extend='both', label=m_diff_base.unit.to_string())
plt.show()

###########################################################################
# Then, we show the running difference map.

plt.figure()
m_diff_running.plot(cmap='Greys_r', norm=norm, title='Running Difference')
plt.colorbar(extend='both', label=m_diff_running.unit.to_string())
plt.show()
