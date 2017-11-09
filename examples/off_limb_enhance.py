"""
===========================
Enhancing Off-limb emission
===========================

This example shows how to enhance emission above the limb.
"""
from __future__ import print_function, division

import numpy as np

import matplotlib.pyplot as plt
import matplotlib.colors as colors

import astropy.units as u

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
# Tell matplotlib to use LaTeX for all the text, and then plot the data and the
# fit.
plt.rc('text', usetex=True)
plt.plot(rsun_array, y, label='data')
label = 'fit=$A\exp(${:.2f}$r)$'.format(params[0])
plt.plot(rsun_array, np.exp(np.poly1d(params)(rsun_array)), label=label)
plt.yscale('log')
plt.ylabel('mean DN')
plt.xlabel('radius r ($R_{\odot}$)')
plt.legend()
plt.show()

###############################################################################
# We now create our scaling array.  At the solar radius, the scale factor is 1.
# Moving away from the disk, the scaling array increases in value.  Finally,
# in order to not affect the emission on the disk, we set the scale factor to
# unity for values of r less than 1.
scale_factor = np.exp((r-1)*-params[0])
scale_factor[r < 1] = 1

###############################################################################
# Let's now plot and compare the results.
scaled_map = sunpy.map.Map(aia.data * scale_factor, aia.meta)
scaled_map.plot_settings['norm'] = colors.Normalize(vmin=10, vmax=10000)

fig = plt.figure(figsize=(12, 5))
ax = fig.add_subplot(121, projection=aia)
aia.plot()
aia.draw_limb()
ax = fig.add_subplot(122, projection=aia)
scaled_map.plot()
scaled_map.draw_limb()
plt.show()
