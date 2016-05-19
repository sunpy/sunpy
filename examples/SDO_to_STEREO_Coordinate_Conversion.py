"""
===================================
AIA to STEREO Coordinate Conversion
===================================

In this example I am going to demonstrate how you can identify a feature on
the surface of the Sun in an AIA image and then convert that point to a point
in a STEREO image.
"""

import matplotlib.pyplot as plt

import astropy.units as u
from astropy.coordinates import SkyCoord

import sunpy.map
import sunpy.coordinates
import sunpy.coordinates.wcs_utils
from sunpy.net import vso


###############################################################################
# The first step is to download some data, we are going to get an image from
# early 2010 when the STEREO spacecraft were roughly 90 deg seperated from the
# Earth.

stereo = (vso.attrs.Source('STEREO_B') &
          vso.attrs.Instrument('EUVI') &
          vso.attrs.Time('2011/1/1', '2011/1/1T00:10:00'))

aia = (vso.attrs.Instrument('AIA') &
       vso.attrs.Sample(24*u.hour) &
       vso.attrs.Time('2011/1/1', '2011/1/2'))

wave = vso.attrs.Wave(30*u.nm, 31*u.nm)


vc = vso.VSOClient()
res = vc.query(wave, aia | stereo)

###############################################################################
# The results from VSO query
print(res)


###############################################################################
# Download the files
files = vc.get(res).wait()


###############################################################################
# Create a dictionary with the two maps, cropped down to full disk.
maps = {m.detector:m.submap((-1100, 1100)*u.arcsec,
                            (-1100, 1100)*u.arcsec) for m in sunpy.map.Map(files)}


###############################################################################
# Plot both maps
fig = plt.figure(figsize=(15,5))
for i, m in enumerate(maps.values()):
    ax = fig.add_subplot(1,2,i+1, projection=m.wcs)
    m.plot(axes=ax)

###############################################################################
# We are now going to pick out a region around the south west corner:

aia_w = 200*u.arcsec
aia_h = 250*u.arcsec
aia_bottom_left = (-800, -300)*u.arcsec


###############################################################################
# Plot a rectangle around the region we want to crop
m = maps['AIA']
fig = plt.figure()
ax = fig.add_subplot(111, projection=m.wcs)
m.plot(axes=ax)
m.draw_rectangle(aia_bottom_left, aia_w, aia_h)


###############################################################################
# Create a submap of this area
subaia = maps['AIA'].submap(u.Quantity((aia_bottom_left[0],
                                        aia_bottom_left[0] + aia_w)),
                            u.Quantity((aia_bottom_left[1],
                                        aia_bottom_left[1] + aia_h)))
subaia.peek()

###############################################################################
# We now want to crop out this same area on the STEREO EUVI image. First, we
# create a `SkyCoord` object with the four corners of the box. When we create
# this object, we use `Map.coordinate_frame` so that the location parameters of
# SDO are correctly set.
hpc_aia = SkyCoord((aia_bottom_left,
                    aia_bottom_left + u.Quantity((aia_w, 0*u.arcsec)),
                    aia_bottom_left + u.Quantity((0*u.arcsec, aia_h)),
                    aia_bottom_left + u.Quantity((aia_w, aia_h))),
                   frame=maps['AIA'].coordinate_frame)

print(hpc_aia)

###############################################################################
# Now we convert these coordinates into Heliographic Stonyhurst coordinates,
# which are on the Sun, with the zero meridian aligned with the Earth.
hgs = hpc_aia.transform_to('heliographic_stonyhurst')
print(hgs)

###############################################################################
# Now we need to provide the position information from the STEREO Imager:
hgs.D0 = maps['EUVI'].dsun
hgs.L0 = maps['EUVI'].heliographic_longitude
hgs.B0 = maps['EUVI'].heliographic_latitude

###############################################################################
# We do this on the Heliographic frame because when in a Heliographic frame
# these parameters have no effect on the frame, but they are used when the
# frame is converted back to Helioprojective. And now we can convert back to
# Helioprojective, but this time from the view point of STEREO B:
hpc_B = hgs.transform_to('helioprojective')
print(hpc_B)

###############################################################################
# Now we can plot this box on both the AIA and EUVI images:
fig = plt.figure(figsize=(15, 5))
for i, (m, coord) in enumerate(zip([maps['EUVI'], maps['AIA']],
                                   [hpc_B, hpc_aia])):
    ax = fig.add_subplot(1, 2, i+1, projection=m.wcs)
    m.plot(axes=ax)

    w = (coord[3].Tx - coord[0].Tx)
    h = (coord[3].Ty - coord[0].Ty)
    m.draw_rectangle(u.Quantity((coord[0].Tx, coord[0].Ty)), w, h,
                     transform=ax.get_transform('world'))

###############################################################################
# We can now zoom in on the region in the EUVI image:
subeuvi = maps['EUVI'].submap(u.Quantity((hpc_B[0].Tx, hpc_B[3].Tx)),
                              u.Quantity((hpc_B[0].Ty, hpc_B[3].Ty)))
subeuvi.peek()

###############################################################################
# Putting them together:
fig = plt.figure(figsize=(15, 5))
for i, m in enumerate((subeuvi, subaia)):
    ax = fig.add_subplot(1, 2, i+1, projection=m.wcs)
    m.plot(axes=ax)

