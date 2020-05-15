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

import astropy.units as u

import sunpy.map
from sunpy.net import Fido
from sunpy.net import attrs as a

######################################################################
# In this example we are going to make a lot of side by side figures, so
# let's change the default figure size.

plt.rcParams['figure.figsize'] = (16, 8)

######################################################################
# We are going to download one AIA and one HMI magnetogram image.

time = (a.Sample(24 * u.hour) &
        a.Time('2010-08-19', '2010-08-19T00:10:00', '2010-08-19') &
        a.vso.Extent(0, 0, 0, 0, "FULLDISK"))
aia = a.Instrument.aia & a.Wavelength(17 * u.nm, 18 * u.nm)
hmi = a.Instrument.hmi & a.Physobs.los_magnetic_field


res = Fido.search(time, aia | hmi)

files = Fido.fetch(res[:, 0])

######################################################################
# We create a map for each image and resample each one just to
# reduce the computation time.

map_aia, map_hmi = [m.resample((1024, 1024)*u.pix) for m in sunpy.map.Map(sorted(files))]
# Why do we have to do this?
map_hmi.plot_settings['cmap'] = "hmimag"
map_hmi.plot_settings['norm'] = plt.Normalize(-2000, 2000)

######################################################################
# Plot both images side by side.

fig = plt.figure()

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
# Construct an output map and set some nice plotting defaults.

out_hmi = sunpy.map.Map(output, map_aia.wcs)
out_hmi.plot_settings['cmap'] = "hmimag"
out_hmi.plot_settings['norm'] = plt.Normalize(-1500, 1500)

fig = plt.figure()
ax1 = fig.add_subplot(1, 2, 1, projection=map_aia)
map_aia.plot(axes=ax1)
ax2 = fig.add_subplot(1, 2, 2, projection=out_hmi)
out_hmi.plot(axes=ax2)


######################################################################
# As both of these images are now on the same pixel grid we can directly plot
# them over one another, by setting the transparency of the HMI plot.

fig = plt.figure(figsize=(20, 15))
ax1 = fig.add_subplot(1, 1, 1, projection=map_aia)
map_aia.plot(axes=ax1)
out_hmi.plot(axes=ax1, alpha=0.5)

plt.show()
