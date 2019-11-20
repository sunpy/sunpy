"""
==========================================
Reprojecting Images to Different Observers
==========================================

This example demonstrates how you can reproject images to the view from
different observers, we use both AIA and STEREO A data to demonstrate this.

You will need `reproject <https://reproject.readthedocs.io/en/stable/>`__ v0.6 or higher installed.
"""
# sphinx_gallery_thumbnail_number = 2

import matplotlib.pyplot as plt
from reproject import reproject_interp

import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.wcs import WCS

import sunpy.map
from sunpy.coordinates import get_body_heliographic_stonyhurst
from sunpy.net import Fido
from sunpy.net import attrs as a

######################################################################
# In this example we are going to make a lot of side by side figures, so
# let's change the default figure size.

plt.rcParams['figure.figsize'] = (16, 8)

######################################################################
# Let’s download an EUV image from both AIA and STEREO A, when their
# separation was around 90 degrees.

stereo = (a.vso.Source('STEREO_A') &
          a.Instrument('EUVI') &
          a.Time('2010-08-19', '2010-08-19T00:10:00'))
aia = (a.Instrument('AIA') &
       a.vso.Sample(24 * u.hour) &
       a.Time('2010-08-19', '2010-08-19T00:10:00'))
wave = a.Wavelength(17 * u.nm, 18 * u.nm)

res = Fido.search(wave, aia | stereo)
files = Fido.fetch(res)

######################################################################
# Create a map for each image.

map_aia, map_stereo = sunpy.map.Map(sorted(files))

fig = plt.figure()
ax1 = fig.add_subplot(1, 2, 1, projection=map_aia)
map_aia.plot(axes=ax1)
ax2 = fig.add_subplot(1, 2, 2, projection=map_stereo)
map_stereo.plot(axes=ax2)

######################################################################
# We now need to construct an output WCS. Because we want to downsample
# the image (for performance) we build a custom header using
# ``sunpy.map.make_fitswcs_header`` but we use a lot of the ``map_aia``
# properties to do it. We use the reference coordinate from the AIA image,
# and make the scale 4x larger to compensate for the fact we are making
# the resolution 4x lower.

out_shape = (1024, 1024)
out_header = sunpy.map.make_fitswcs_header(
    out_shape,
    map_aia.reference_coordinate,
    scale=u.Quantity(map_aia.scale)*4,
    instrument="EUVI",
    observatory="AIA Observer",
    wavelength=map_stereo.wavelength
)

######################################################################
# Next we construct an `~astropy.wcs.WCS` object from the header.
# Currently `~astropy.wcs.WCS` does not understand the observer
# position, so we manually set that.

out_wcs = WCS(out_header)
out_wcs.heliographic_observer = map_aia.reference_coordinate.observer

######################################################################
# We can now reproject the STEREO map to this output `~astropy.wcs.WCS`.
# Here we are using the fastest but least accurate method of reprojection,
# `reproject.reproject_interp`, a more accurate but slower method is
# `reproject.reproject_adaptive`.

output, footprint = reproject_interp(map_stereo, out_wcs, out_shape)

######################################################################
# We can now plot the STEREO image as seen from the position of SDO, next
# to the AIA image.

outmap = sunpy.map.Map(output, out_header)
outmap.plot_settings = map_stereo.plot_settings

fig = plt.figure()
ax1 = fig.add_subplot(1, 2, 1, projection=map_aia)
map_aia.plot(axes=ax1)
ax2 = fig.add_subplot(1, 2, 2, projection=outmap)
outmap.plot(axes=ax2)

######################################################################
# AIA as Seen from Mars
# =====================
#
# We can also change the observer of the AIA image to any observer
# coordinate. SunPy provides a function which can get the observer
# coordinate for any known body. In this example we are going to use Mars.

mars = get_body_heliographic_stonyhurst('mars', map_aia.date)

######################################################################
# We now generate a target WCS, to do this we need to generate a reference
# coordinate, which is similar to the aia frame, but scaled down by 4x and
# with the observer at Mars. To do this we generate a new reference
# coordinate.

mars_ref_coord = SkyCoord(map_aia.reference_coordinate.Tx,
                          map_aia.reference_coordinate.Ty,
                          obstime=map_aia.reference_coordinate.obstime,
                          observer=mars,
                          frame="helioprojective")

######################################################################
# then a header

out_shape = (1024, 1024)
mars_header = sunpy.map.make_fitswcs_header(
    out_shape,
    mars_ref_coord,
    scale=u.Quantity(map_aia.scale)*4,
    rotation_matrix=map_aia.rotation_matrix,
    instrument="AIA",
    wavelength=map_aia.wavelength
)

######################################################################
# Once again we need to generate a `~astropy.wcs.WCS` and then manually
# set the observer location.

mars_wcs = WCS(out_header)
mars_wcs.heliographic_observer = mars

output, footprint = reproject_interp(map_aia, mars_wcs, out_shape)

######################################################################
# We generate the output map and plot it next to the original image.

outmap = sunpy.map.Map((output, mars_header))
outmap.plot_settings = map_aia.plot_settings

fig = plt.figure()

ax1 = fig.add_subplot(1, 2, 1, projection=map_aia)
map_aia.plot(axes=ax1)
outmap.draw_grid(color='w')

ax2 = fig.add_subplot(1, 2, 2, projection=outmap)
outmap.plot(axes=ax2)
outmap.draw_grid(color='w')

plt.show()
