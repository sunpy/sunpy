"""
===========================
Plotting a difference image
===========================

How to compute and plot a difference image.
The example uses `sunpy.map.MapSequence` to compute a difference image and then plot it.
This basic method works for base difference or running difference. Just change whether
you're subtracting the previous image or the first image in a sequence.
"""
# sphinx_gallery_thumbnail_number = 2
import matplotlib.colors as colors
import matplotlib.pyplot as plt

import astropy.units as u

import sunpy.data.sample
import sunpy.map

###########################################################################
# First load a couple of images taken from the sample dataset. These are
# two cutouts taken during a flare.
maps = sunpy.map.Map([sunpy.data.sample.AIA_193_CUTOUT03_IMAGE,
                      sunpy.data.sample.AIA_193_CUTOUT04_IMAGE],
                     sequence=True)

###########################################################################
# Now we'll do a standard plot of the second image just to see it.

plt.figure()
ax = plt.subplot(projection=maps[1])
maps[1].plot(clip_interval=(0.5, 99.9)*u.percent)

###########################################################################
# And now we can do take the actual difference.

diff = maps[1].data - maps[0].data

###########################################################################
# But we have to decide what to do with the metadata. For example, what
# time does this difference image correspond to? The time of the first or
# second image? The mean time? You'll have to decide what makes most sense
# for your application. Here we'll just use the metadata from the second
# image. Then we can store the difference and header back in a Map.

meta = maps[1].meta
diff_map = sunpy.map.Map(diff, meta)

###########################################################################
# Finally, we'll plot it. We'll apply a colormap and re-normalize the
# intensity so that it shows up well.

plt.figure()
plt.subplot(projection=diff_map)
diff_map.plot(cmap='Greys_r',
              norm=colors.Normalize(vmin=-200, vmax=200))
plt.colorbar(extend='both', label=maps[1].unit.to_string())
plt.show()
