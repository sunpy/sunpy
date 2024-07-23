"""
=========================================
Tracking and Co-aligning an Active Region
=========================================

This example demonstrates how to track a particular region as a function of time,
create a cutout around that region, and align it at each time step to get an aligned datacube.
"""

import matplotlib
import matplotlib.pyplot as plt

matplotlib.use('Agg')  # Use a non-interactive backend suitable for script and doc builds

import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.visualization import ImageNormalize, SqrtStretch

import sunpy.map
from sunpy.coordinates import propagate_with_solar_surface
from sunpy.net import Fido
from sunpy.net import attrs as a

###############################################################################
# As we require a series of images, we will need to download them using `sunpy.net.Fido`.
# For this example, we will acquire 4 images at 6-hour intervals.

query = Fido.search(a.Time('2012/3/4', '2012/3/5'),
                    a.Instrument.aia,
                    a.Wavelength(171*u.angstrom),
                    a.Sample(360*u.minute))
print(query)
files = Fido.fetch(query)

###############################################################################
# Now we have a set of images, we can create a sequence of maps.

aia_sequence = sunpy.map.Map(files, sequence=True)

fig = plt.figure()
ax = fig.add_subplot(projection=aia_sequence.maps[0])
ani = aia_sequence.plot(axes=ax, norm=ImageNormalize(vmin=0, vmax=5e3, stretch=SqrtStretch()))

plt.show()

###############################################################################
# Now, let us crop into an interesting region.
corner = SkyCoord(Tx=-375*u.arcsec, Ty=0*u.arcsec, frame=aia_sequence[0].coordinate_frame)

# Create a cutout out from the bottom left corner.
cutout_map = aia_sequence[0].submap(corner, width=500*u.arcsec, height=500*u.arcsec)

fig = plt.figure()
ax = fig.add_subplot(projection=cutout_map)

cutout_map.plot(axes=ax)

###############################################################################
# Track and co-align the region across the sequence of maps using solar rotation.
fig = plt.figure(figsize=(24, 8))
for i, m in enumerate(aia_sequence):
    ax = fig.add_subplot(1, len(aia_sequence), i+1, projection=m)
    m.plot(axes=ax)
    ax.set_xlabel(' ')
    ax.set_ylabel(' ')
    plt.subplots_adjust(wspace=0.3)
    with propagate_with_solar_surface():
        blc = cutout_map.bottom_left_coord.transform_to(m.coordinate_frame)
        trc = cutout_map.top_right_coord.transform_to(m.coordinate_frame)
        m.draw_quadrangle(blc, top_right=trc)

###############################################################################
# Aligns all images to the World Coordinate System (WCS) of the cutout.
# `propagate_with_solar_surface` is a context manager that adjusts for solar rotation,
# ensuring that regions on the Sun's surface remain in the correct position as the Sun rotates.

with propagate_with_solar_surface():
    aia_sequence_aligned = sunpy.map.Map([m.reproject_to(cutout_map.wcs) for m in aia_sequence], sequence=True)

# `reproject_to` reprojects each map in the sequence to the WCS of the cutout map,
# aligning all images to the same reference frame.

# Plot the aligned sequence of maps.
fig = plt.figure(figsize=(24, 8))
for i, m in enumerate(aia_sequence_aligned):
    ax = fig.add_subplot(1, len(aia_sequence_aligned), i+1, projection=m)
    m.plot(axes=ax, cmap='sdoaia171', title=aia_sequence[i].date)
    ax.set_xlabel(' ')
    ax.set_ylabel(' ')
plt.subplots_adjust(wspace=0.3)
# This code aligns a sequence of solar images to a common reference frame, taking into account solar rotation,
# and then plots them side by side in a single figure. Each subplot shows one map from the aligned sequence.
