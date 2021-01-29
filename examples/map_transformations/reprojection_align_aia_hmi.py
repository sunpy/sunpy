"""
========================================
Aligning AIA and HMI Data with Reproject
========================================

A common case when using data from multiple sources,
is aligning so they have the same reference frame.
It is also possible to use `reproject` to align data, by reprojecting one image
to the WCS of another. This is a very generic way of aligning data, and can be
very accurate.

You will need `reproject <https://reproject.readthedocs.io/en/stable/>`__ v0.6 or higher installed.
"""
import matplotlib.pyplot as plt
from reproject import reproject_interp

import sunpy.data.sample
import sunpy.map

######################################################################
# We use the AIA image and HMI image from the sample data.  For the
# HMI map, we use the special HMI color map, which expects the plotted
# range to be -1500 to 1500.

map_aia = sunpy.map.Map(sunpy.data.sample.AIA_171_IMAGE)
map_hmi = sunpy.map.Map(sunpy.data.sample.HMI_LOS_IMAGE)
map_hmi.plot_settings['cmap'] = "hmimag"
map_hmi.plot_settings['norm'] = plt.Normalize(-1500, 1500)

######################################################################
# Plot both images side by side.

fig = plt.figure(figsize=(12, 5))
ax1 = fig.add_subplot(1, 2, 1, projection=map_aia)
map_aia.plot(axes=ax1)
ax2 = fig.add_subplot(1, 2, 2, projection=map_hmi)
map_hmi.plot(axes=ax2)

######################################################################
# We can now reproject the HMI image to the WCS of the AIA image. We are using
# the fast `~reproject.reproject_interp`, however the slower but most accurate
# `~reproject.reproject_exact` would also work well here. The
# `~reproject.reproject_exact` function only works when reprojecting between
# two WCSes with the same observer, which makes it well suited to aligning
# data.

output, footprint = reproject_interp(map_hmi, map_aia.wcs, map_aia.data.shape)

######################################################################
# Construct an output map and use the same HMI plot settings.

out_hmi = sunpy.map.Map(output, map_aia.wcs)
out_hmi.plot_settings = map_hmi.plot_settings

fig = plt.figure(figsize=(12, 5))
ax1 = fig.add_subplot(1, 2, 1, projection=map_aia)
map_aia.plot(axes=ax1)
ax2 = fig.add_subplot(1, 2, 2, projection=out_hmi)
out_hmi.plot(axes=ax2, title='Reprojected HMI image')

######################################################################
# As both of these images are now on the same pixel grid we can directly plot
# them over one another, by setting the transparency of the HMI plot.

fig = plt.figure()
ax1 = fig.add_subplot(1, 1, 1, projection=map_aia)
map_aia.plot(axes=ax1)
out_hmi.plot(axes=ax1, alpha=0.5)
plt.title('HMI overlaid on AIA')

plt.show()
