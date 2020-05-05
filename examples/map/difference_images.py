"""
===========================
Plotting a difference image
===========================

How to compute and plot a difference image.
The example uses `sunpy.map.MapSequence` to compute a difference image and then plot it.
This basic method works for base difference or running difference. Just change whether
you're subtracting the previous image or the first image in a sequence.
"""
###########################################################################

import matplotlib.colors as colors
import matplotlib.pyplot as plt

import astropy.units as u

import sunpy.map
from sunpy.net import Fido
from sunpy.net import attrs as a

###########################################################################
# First we'll download a couple images and store them in a
# `sunpy.map.MapSequence`. This could be from any channel of any imager.
# Here, we use SDO/AIA 304 Å.

instrument = a.Instrument.aia
wave = a.Wavelength(30 * u.nm, 31 * u.nm)
result = Fido.search(a.Time('2015-06-18T00:00:00', '2015-06-18T00:00:10') |
                     a.Time('2015-06-18T01:03:30', '2015-06-18T01:03:35'),
                     instrument,
                     wave)
downloaded_files = Fido.fetch(result)
maps = sunpy.map.Map(downloaded_files, sequence=True)

###########################################################################
# Now we'll do a standard plot of the second image just to see it.

plt.figure()
ax = plt.subplot(projection=maps[1])
maps[1].plot()

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
# Finally, we'll plot it. We'll apply a colormap and renormalize the
# intensity so that it shows up well.

plt.figure()
ax_diff = plt.subplot(projection=diff_map)
diff_map.plot(cmap='Greys_r',
              norm=colors.Normalize(vmin=-50, vmax=50))
plt.show()

# sphinx_gallery_thumbnail_number = 2
