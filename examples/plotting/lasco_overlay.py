"""
==============================================
Overlay an AIA image on a LASCO C2 coronagraph
==============================================

This example shows the steps needed to overlay the disc and off-limb
components of an AIA image within the masked occulter of a LASCO C2 image.

It also demonstrates the usage of ``hypy`` to download data from the Helioviewer Project.
"""
from datetime import datetime

import hvpy
import matplotlib.pyplot as plt
from hvpy.datasource import DataSource

import astropy.units as u
from astropy.coordinates import SkyCoord

import sunpy.data.sample
from sunpy.coordinates import Helioprojective
from sunpy.map import Map

###############################################################################
# First, we will acquire a calibrated LASCO C2 image from Helioviewer and
# create a map. hvpy uses the standard datetime instead of astropy.time.

lasco_jp2_file = hvpy.save_file(hvpy.getJP2Image(datetime(2011, 6, 7, 6, 34),
                                                 DataSource.LASCO_C2.value),
                                filename="lasco.jp2", overwrite=True)
lasco_map = Map(lasco_jp2_file)
aia_map = Map(sunpy.data.sample.AIA_171_IMAGE)

###############################################################################
# In order to plot off-limb features of the AIA image, we need to reproject AIA
# using `~sunpy.coordinates.Helioprojective.assume_spherical_screen`.

projected_coord = SkyCoord(0*u.arcsec, 0*u.arcsec,
                           obstime=lasco_map.observer_coordinate.obstime,
                           frame='helioprojective',
                           observer=lasco_map.observer_coordinate,
                           rsun=aia_map.coordinate_frame.rsun)
projected_header = sunpy.map.make_fitswcs_header(aia_map.data.shape,
                                                 projected_coord,
                                                 scale=u.Quantity(aia_map.scale),
                                                 instrument=aia_map.instrument,
                                                 wavelength=aia_map.wavelength)
with Helioprojective.assume_spherical_screen(aia_map.observer_coordinate):
    aia_reprojected = aia_map.reproject_to(projected_header)

###############################################################################
# Finally, we plot the images by layering the AIA image on top of the LASCO C2
# image.

fig = plt.figure()
ax = fig.add_subplot(projection=lasco_map)
lasco_map.plot(axes=ax)
aia_reprojected.plot(axes=ax, clip_interval=(1, 99.9)*u.percent, autoalign=True)
ax.set_title("AIA and LASCO C2 Overlay")

plt.show()
