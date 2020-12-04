"""
===================================
AIA to STEREO coordinate conversion
===================================

How to convert a point of a source on an AIA image
to a position on a STEREO image.
"""
import matplotlib.pyplot as plt

import astropy.units as u
from astropy.coordinates import SkyCoord

import sunpy.coordinates
import sunpy.map
from sunpy.net import Fido
from sunpy.net import attrs as a

###############################################################################
# The first step is to download some data, we are going to get an image from
# early 2011 when the STEREO spacecraft were roughly 90 deg seperated from the
# Earth.
stereo = (a.Source('STEREO_B') &
          a.Instrument("EUVI") &
          a.Time('2011-01-01', '2011-01-01T00:10:00'))

aia = (a.Instrument.aia &
       a.Sample(24 * u.hour) &
       a.Time('2011-01-01', '2011-01-02'))

wave = a.Wavelength(30 * u.nm, 31 * u.nm)


res = Fido.search(wave, aia | stereo)

###############################################################################
# Download the files:
files = Fido.fetch(res)
print(files)

###############################################################################
# Create a dictionary with the two maps, cropped down to full disk.
maps = {m.detector: m.submap(SkyCoord([-1100, 1100]*u.arcsec,
                                      [-1100, 1100]*u.arcsec,
                                      frame=m.coordinate_frame))
        for m in sunpy.map.Map(files)}

###############################################################################
# Plot both maps
fig = plt.figure(figsize=(10, 4))
for i, m in enumerate(maps.values()):
    ax = fig.add_subplot(1, 2, i+1, projection=m)
    m.plot(axes=ax)

###############################################################################
# We are now going to pick out a region around the south west corner:
aia_bottom_left = SkyCoord(-800 * u.arcsec,
                           -300 * u.arcsec,
                           frame=maps['AIA'].coordinate_frame)
aia_top_right = SkyCoord(-600 * u.arcsec,
                         -50 * u.arcsec,
                         frame=maps['AIA'].coordinate_frame)

###############################################################################
# Plot a rectangle around the region we want to crop
m = maps['AIA']
fig = plt.figure()
ax = fig.add_subplot(111, projection=m)
m.plot(axes=ax)
m.draw_rectangle(aia_bottom_left, top_right=aia_top_right)


###############################################################################
# Create a submap of this area
subaia = maps['AIA'].submap(aia_bottom_left, top_right=aia_top_right)
fig = plt.figure()
subaia.plot()

###############################################################################
# We now want to crop out this same area on the STEREO EUVI image. First, we
# create a `SkyCoord` object with the four corners of the box. When we create
# this object, we use `Map.coordinate_frame` so that the location parameters of
# SDO are correctly set.
corners = ([aia_bottom_left.Tx, aia_bottom_left.Ty],
           [aia_top_right.Tx, aia_bottom_left.Ty],
           [aia_bottom_left.Tx, aia_top_right.Ty],
           [aia_top_right.Tx, aia_top_right.Ty])

hpc_aia = SkyCoord(corners, frame=maps['AIA'].coordinate_frame)

print(hpc_aia)

###############################################################################
# We can now transform to from the AIA frame to the EUVI frame.
# This transformation first transforms to Heliographic Stonyhurst coordinates
# and then into the EUVI frame.
hpc_B = hpc_aia.transform_to(maps['EUVI'].coordinate_frame)
print(hpc_B)

###############################################################################
# Now we can plot this box on both the AIA and EUVI images:
fig = plt.figure(figsize=(10, 4))
for i, (m, coord) in enumerate(zip([maps['EUVI'], maps['AIA']],
                                   [hpc_B, hpc_aia])):
    ax = fig.add_subplot(1, 2, i+1, projection=m)
    m.plot(axes=ax)

    # coord[3] is the top-right corner coord[0] is the bottom-left corner.
    m.draw_rectangle(coord[0], top_right=coord[3],
                     transform=ax.get_transform('world'))

###############################################################################
# We can now zoom in on the region in the EUVI image:
subeuvi = maps['EUVI'].submap(hpc_B[0], top_right=hpc_B[3])
fig = plt.figure()
plt.subplot(projection=subeuvi)
subeuvi.plot()

###############################################################################
# Putting them together:
fig = plt.figure(figsize=(15, 5))
for i, m in enumerate((subeuvi, subaia)):
    ax = fig.add_subplot(1, 2, i+1, projection=m)
    m.plot(axes=ax)
plt.show()
