"""
===========================
Enhancing Off-limb emission
===========================

This example shows how to enhance emission above the limb.
"""
import numpy as np
import astropy.units as u
import matplotlib.pyplot as plt
import matplotlib.colors as colors

from sunpy.data.sample import AIA_171_IMAGE
from sunpy.map import Map

###############################################################################
# We first create the Map using the sample data.
aia = Map(AIA_171_IMAGE)

###############################################################################
# Next we build two arrays which include all of the x and y pixel indices.
# We must not forget to add the correct units because we will next pass this
# into a function which requires them.
x, y = np.meshgrid(*[np.arange(v.value) for v in aia.dimensions]) * u.pix

###############################################################################
# Now we can convert this to helioprojective coordinates and create a new
# array which contains the normalized radial position for each pixel
hpc_coords = aia.pixel_to_data(x,y)
r = np.sqrt(hpc_coords.Tx ** 2 + hpc_coords.Ty ** 2) / aia.rsun_obs

###############################################################################
# Let's check how emission above the limb depends on distance
rsun_step_size = 0.01
rsun_array = np.arange(1, r.max(), rsun_step_size)
y = np.array([aia.data[(r > this_r) * (r < this_r + rsun_step_size)].mean() for this_r in rsun_array])

###############################################################################
# Next let's plot it along with a fit to the data. We perform the fit in
# log-linear space.
params = np.polyfit(rsun_array[rsun_array<1.5], np.log10(y[rsun_array<1.5]), 1)

plt.plot(rsun_array, y, label='data')
plt.plot(rsun_array, 10**np.poly1d(params)(rsun_array), label=r'fit=Ax^{:.2f}'.format(params[0]))
plt.yscale('log')
plt.ylabel('mean DN')
plt.xlabel('Rsun')
plt.legend()
plt.savefig("gallery_offlimb_dependence.png")
#plt.show()
###############################################################################
# We now create our scaling array which we will multiply our original data by 5.
# In order to not affect the emission on the disk, we set the scale factor to
# unity for values of r below 1.
scale_exponent = 5
scale_factor =  r ** scale_exponent
scale_factor[r < 1] = 1

###############################################################################
# Let's now plot and compare the results.
scaled_map = Map(aia.data * scale_factor, aia.meta)
norm = colors.Normalize(vmin=10, vmax=10000)
plt.figure(figsize=(12,5))
plt.subplot(1,2,1)
aia.plot()
plt.subplot(1,2,2)
scaled_map.plot()
plt.savefig('gallery_offlimb_enhance.png')
#plt.show()
