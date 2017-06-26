"""
=============
Map Histogram
=============

How to inspect the histogram of the data of a map.
"""
import numpy as np
import matplotlib.pyplot as plt
import astropy.units as u
from astropy.coordinates import SkyCoord
from sunpy.data.sample import AIA_171_IMAGE
import sunpy.map

###############################################################################
# We first create the Map using the sample data and we will create a submap
# of a quiet region.
aia = sunpy.map.Map(AIA_171_IMAGE)
bl = SkyCoord(-400 * u.arcsec, 0 * u.arcsec, frame=aia.coordinate_frame)
tr = SkyCoord(0 * u.arcsec, 400 * u.arcsec, frame=aia.coordinate_frame)
aia_smap = aia.submap(bl, tr)

###############################################################################
# We now create a histogram of the data in this region.
dmin = aia_smap.min()
dmax = aia_smap.max()
num_bins = 50
hist, bins = np.histogram(aia_smap.data, bins=np.linspace(dmin, dmax, num_bins))
width = 0.7 * (bins[1] - bins[0])
x = (bins[:-1] + bins[1:]) / 2

###############################################################################
# Let's plot the histogram as well as some standard values such as mean
# upper, and lower value and the one-sigma range.
plt.figure()
plt.bar(x, hist, align='center', width=width, label='Histogram')
plt.xlabel('Intensity')
plt.axvline(dmin, label='Data min={:.2f}'.format(dmin), color='black')
plt.axvline(dmax, label='Data max={:.2f}'.format(dmax), color='black')
plt.axvline(aia_smap.data.mean(),
            label='mean={:.2f}'.format(aia_smap.data.mean()), color='green')
one_sigma = np.array([aia_smap.data.mean() - aia_smap.data.std(),
                      aia_smap.data.mean() + aia_smap.data.std()])
plt.axvspan(one_sigma[0], one_sigma[1], alpha=0.3, color='green',
            label='mean +/- std = [{0:.2f}, {1:.2f}]'.format(
            one_sigma[0], one_sigma[1]))
plt.axvline(one_sigma[0], color='green')
plt.axvline(one_sigma[1], color='red')
plt.legend()
plt.show()

###############################################################################
# Finally let's overplot what the one-sigma range means on the map
fig = plt.figure()
fig.add_subplot(projection=aia_smap)
aia_smap.plot()
levels = one_sigma / dmax * u.percent * 100
aia_smap.draw_contours(levels=levels, colors=['red', 'green'])
plt.show()
