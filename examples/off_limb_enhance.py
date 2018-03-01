"""
===========================
Enhancing Off-limb emission
===========================

This example shows how to enhance emission above the limb.
"""
from __future__ import print_function, division

import numpy as np

import matplotlib.pyplot as plt

import astropy.units as u
from astropy.visualization.mpl_normalize import ImageNormalize

import sunpy.map
from sunpy.data.sample import AIA_171_IMAGE

###############################################################################
# We first create the Map using the sample data.

aia = sunpy.map.Map(AIA_171_IMAGE)

###############################################################################
# Next we build two arrays which include all of the x and y pixel indices.
# We must not forget to add the correct units because we will next pass this
# into a function which requires them.

x, y = np.meshgrid(*[np.arange(v.value) for v in aia.dimensions]) * u.pix

###############################################################################
# Now we can convert this to helioprojective coordinates and create a new
# array which contains the normalized radial position for each pixel

hpc_coords = aia.pixel_to_world(x, y)
r = np.sqrt(hpc_coords.Tx ** 2 + hpc_coords.Ty ** 2) / aia.rsun_obs

###############################################################################
# Let's check how emission above the limb depends on distance

rsun_step_size = 0.01
rsun_array = np.arange(1, r.max(), rsun_step_size)
y = np.array([aia.data[(r > this_r) * (r < this_r + rsun_step_size)].mean()
              for this_r in rsun_array])

###############################################################################
# Next let's plot it along with a fit to the data. We perform the fit in
# linear-log space.  We fit the logarithm of the intensity since the intensity
# drops of very quickly as a function of distance from the limb.

params = np.polyfit(rsun_array[rsun_array < 1.5],
                    np.log(y[rsun_array < 1.5]), 1)

###############################################################################
# Tell matplotlib to use LaTeX for all the text, make the fontsize bigger, and
# then plot the data and the fit.

fontsize = 14
plt.plot(rsun_array, y, label='data')
best_fit = np.exp(np.poly1d(params)(rsun_array))
label = r'best fit: {:.2f}$e^{{{:.2f}r}}$'.format(best_fit[0], params[0])
plt.plot(rsun_array, best_fit, label=label)
plt.yscale('log')
plt.ylabel(r'mean DN', fontsize=fontsize)
plt.xlabel(r'radius r ($R_{\odot}$)', fontsize=fontsize)
plt.xticks(fontsize=fontsize)
plt.yticks(fontsize=fontsize)
plt.title(r'observed off limb mean DN and best fit', fontsize=fontsize)
plt.legend(fontsize=fontsize)
plt.tight_layout()
plt.show()

###############################################################################
# We now create our scaling array.  At the solar radius, the scale factor is 1.
# Moving away from the disk, the scaling array increases in value.  Finally,
# in order to not affect the emission on the disk, we set the scale factor to
# unity for values of r less than 1.

scale_factor = np.exp((r-1)*-params[0])
scale_factor[r < 1] = 1

###############################################################################
# Let's now plot and compare the results.  The scaled map uses the same image
# stretching function as the original image (set by the keyword 'stretch')
# clipped to the same range (set by the keywords 'vmin' and 'vmax').

scaled_map = sunpy.map.Map(aia.data * scale_factor, aia.meta)
scaled_map.plot_settings['norm'] = ImageNormalize(stretch=aia.plot_settings['norm'].stretch,
                                                  vmin=aia.data.min(), vmax=aia.data.max())

fig = plt.figure(figsize=(12, 5))
ax = fig.add_subplot(121, projection=aia)
aia.plot()
aia.draw_limb()
ax = fig.add_subplot(122, projection=aia)
scaled_map.plot()
scaled_map.draw_limb()
plt.show()
