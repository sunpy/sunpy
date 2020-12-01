"""
===============================================
Requesting cutouts of AIA images from the JSOC
===============================================

This example shows how to request a cutout of a series of
AIA images from the JSOC and animate the resulting sequence.
"""
import os

import matplotlib.pyplot as plt

import astropy.time
import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.visualization import ImageNormalize, SqrtStretch

import sunpy.map
from sunpy.net import Fido
from sunpy.net import attrs as a
from sunpy.net.jsoc import JSOCClient

#####################################################
# First, query a full frame AIA image.
t0 = astropy.time.Time('2012-09-24T14:56:03', scale='utc', format='isot')
q = Fido.search(
    a.Instrument.aia,
    a.Physobs.intensity,
    a.Wavelength(171*u.angstrom),
    a.Time(t0, t0 + 13*u.s),
)
m = sunpy.map.Map(Fido.fetch(q))

#####################################################
# Next, we will create a submap from this image. We will
# crop the field of view to active region NOAA 11575.
m_cutout = m.submap(
    SkyCoord(-500*u.arcsec, -275*u.arcsec, frame=m.coordinate_frame),
    top_right=SkyCoord(150*u.arcsec, 375*u.arcsec, frame=m.coordinate_frame),
)
m_cutout.peek()

#####################################################
# We want to watch the evolution of this active region
# in time, but we do not want to download the full frame
# image at each timestep. Instead, we will use our submap
# to create a cutout request from the JSOC.
#
# First, construct the cutout from the submap
# above using the `~sunpy.net.jsoc.attrs.Cutout` attribute.
cutout = a.jsoc.Cutout(
    m_cutout.bottom_left_coord,
    top_right=m_cutout.top_right_coord,
    tracking=True
)

#####################################################
# Now we are ready to construct the query. Note that all of this is
# the same for a full-frame image except for the
# cutout component.
c = JSOCClient()
q = c.search(
    a.Time(m_cutout.date - 6*u.h, m_cutout.date + 6*u.h),
    a.Wavelength(m_cutout.wavelength),
    a.Sample(12*u.min),
    a.jsoc.Series.aia_lev1_euv_12s,
    a.jsoc.Notify(os.environ["JSOC_EMAIL"]),  # Put your email here
    a.jsoc.Segment.image,
    cutout,
)

#####################################################
# Submit the export request and download the data. This
# may take some time as the data has to be processed and
# downloaded.
req = c.request_data(q, method='url', protocol='fits')
while not req.has_succeeded():
    continue
files = c.get_request(req, max_conn=2)
files.sort()

#####################################################
# Now that we've downloaded the files, we can create
# a `~sunpy.map.MapSequence` from them.
m_seq = sunpy.map.Map(files, sequence=True)

#####################################################
# Finally, we can construct an animation in time from
# our stack of cutouts and interactively flip through
# each image in our sequence.
for m in m_seq:
    m.plot_settings['norm'] = ImageNormalize(vmin=0, vmax=5e3, stretch=SqrtStretch())
m_seq.peek()
plt.show()
