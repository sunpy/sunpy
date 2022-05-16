"""
===============================================
Requesting cutouts of AIA images from the JSOC
===============================================

This example shows how to request a cutout of a series of
AIA images from the JSOC and animate the resulting sequence.
"""
# sphinx_gallery_thumbnail_number = 2
import os

import matplotlib.pyplot as plt

import astropy.time
import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.visualization import ImageNormalize, SqrtStretch

import sunpy.map
from sunpy.net import Fido
from sunpy.net import attrs as a

#####################################################
# First, query a full frame AIA image.

t0 = astropy.time.Time('2012-09-24T14:56:03', scale='utc', format='isot')
query = Fido.search(
    a.Instrument.aia,
    a.Physobs.intensity,
    a.Wavelength(171*u.angstrom),
    a.Time(t0, t0 + 13*u.s),
)
files = Fido.fetch(query)
amap = sunpy.map.Map(files)

#####################################################
# Next, we will create a submap from this image. We will
# crop the field of view to active region NOAA 11575.

cutouts = amap.submap(
    SkyCoord(-500*u.arcsec, -275*u.arcsec, frame=amap.coordinate_frame),
    top_right=SkyCoord(150*u.arcsec, 375*u.arcsec, frame=amap.coordinate_frame),
)
plt.figure()
cutouts.plot()

plt.show()

#####################################################
# We want to watch the evolution of this active region
# in time, but we do not want to download the full frame
# image at each timestep. Instead, we will use our submap
# to create a cutout request from the JSOC.
#
# First, construct the cutout from the submap
# above using the `~sunpy.net.jsoc.attrs.Cutout` attribute.

cutout = a.jsoc.Cutout(
    cutouts.bottom_left_coord,
    top_right=cutouts.top_right_coord,
    tracking=True
)

#####################################################
# Exporting data from the JSOC requires registering your
# email first. Please replace this with your email
# address once you have registered.
# See `this page <http://jsoc.stanford.edu/ajax/register_email.html>`_
# for more details.

jsoc_email = os.environ["JSOC_EMAIL"]

#####################################################
# Now we are ready to construct the query. Note that all of this is
# the same for a full-frame image except for the
# cutout component. We will download images from a 12 hour interval
# centered on the time of the above cutout.
# We request one image every 2 hours.

query = Fido.search(
    a.Time(cutouts.date - 6*u.h, cutouts.date + 6*u.h),
    a.Wavelength(cutouts.wavelength),
    a.Sample(2*u.h),
    a.jsoc.Series.aia_lev1_euv_12s,
    a.jsoc.Notify(jsoc_email),
    a.jsoc.Segment.image,
    cutout,
)

#####################################################
# Submit the export request and download the data.

files = Fido.fetch(query)
files.sort()

#####################################################
# Now that we've downloaded the files, we can create
# a `~sunpy.map.MapSequence` from them.

sequence = sunpy.map.Map(files, sequence=True)

#####################################################
# Finally, we can construct an animation in time from
# our stack of cutouts and interactively flip through
# each image in our sequence. We first adjust the plot
# settings on each image to ensure the colorbar is the
# same at each time step.

for each_map in sequence:
    each_map.plot_settings['norm'] = ImageNormalize(vmin=0, vmax=5e3, stretch=SqrtStretch())
plt.figure()
ani = sequence.plot()

plt.show()
