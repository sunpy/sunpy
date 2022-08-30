"""
===============================================
Requesting cutouts of AIA images from the JSOC
===============================================

This example shows how to request a cutout of a series of
AIA images from the JSOC.
"""
# sphinx_gallery_thumbnail_number = 2
import os

import astropy.time
import astropy.units as u
from astropy.coordinates import SkyCoord

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
# Next, we will use this map to definte the top right and bottom
# left coordinates we want for the submap.
bottom_left = SkyCoord(-500*u.arcsec, -275*u.arcsec, frame=amap.coordinate_frame)
top_right = SkyCoord(150*u.arcsec, 375*u.arcsec, frame=amap.coordinate_frame)

#####################################################
# Now construct the cutout from the coordinates above
# above using the `~sunpy.net.jsoc.attrs.Cutout` attribute.

cutout = a.jsoc.Cutout(bottom_left, top_right=top_right, tracking=True)

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
    a.Time(amap.date - 6*u.h, amap.date + 6*u.h),
    a.Wavelength(amap.wavelength),
    a.Sample(2*u.h),
    a.jsoc.Series.aia_lev1_euv_12s,
    a.jsoc.Notify(jsoc_email),
    a.jsoc.Segment.image,
    cutout,
)
print(query)

#####################################################
# To download use the following code:
# ``files = Fido.fetch(query)``
