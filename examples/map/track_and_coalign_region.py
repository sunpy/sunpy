"""
=========================================
Tracking and Co-aligning an Active Region
=========================================

This example demonstrates how to track a particular region as a function of time,
create a cutout around that region, and align it at each time step to get an aligned datacube.
"""

import matplotlib.pyplot as plt

import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.visualization import ImageNormalize, SqrtStretch

import sunpy.map
from sunpy.coordinates import propagate_with_solar_surface
from sunpy.net import Fido
from sunpy.net import attrs as a

###############################################################################
# First, let's download a series of images in time using `~sunpy.net.Fido`.
# In this example, we will download a series of AIA 171 Å images observed over the course
# of one day at a cadence of 1 image every 6 hours.

query = Fido.search(a.Time('2012/3/4', '2012/3/5'),
                    a.Instrument.aia,
                    a.Wavelength(171*u.angstrom),
                    a.Sample(6*u.h))
print(query)
files = Fido.fetch(query)

###############################################################################
# Now that we have a set of images in time, we can create a `~sunpy.map.MapSequence` to hold all of them
# and animate that sequence in time.

aia_sequence = sunpy.map.Map(files, sequence=True)

fig = plt.figure()
ax = fig.add_subplot(projection=aia_sequence[0])
ani = aia_sequence.plot(axes=ax, norm=ImageNormalize(vmin=0, vmax=5e3, stretch=SqrtStretch()))

###############################################################################
# Next, let's crop one of the maps in our sequence to the active region of interest.

corner = SkyCoord(Tx=280*u.arcsec, Ty=0*u.arcsec, frame=aia_sequence[0].coordinate_frame)
cutout_map = aia_sequence[0].submap(corner, width=500*u.arcsec, height=500*u.arcsec)

fig = plt.figure()
ax = fig.add_subplot(projection=cutout_map)
cutout_map.plot(axes=ax)

###############################################################################
# Now that we have our cutout around our active region, we can reproject each map in our sequence to
# the WCS of that cutout. Additionally, we will use the `~sunpy.coordinates.propagate_with_solar_surface`
# context manager to adjust the field of view of the cutout with the rotation of the solar surface.

with propagate_with_solar_surface():
    aia_sequence_aligned = sunpy.map.Map([m.reproject_to(cutout_map.wcs) for m in aia_sequence], sequence=True)

###############################################################################
# Finally, we can animate our sequence of reprojected cutouts to confirm that we've tracked our active
# region of interest as a function of time across the solar disk.

fig = plt.figure()
ax = fig.add_subplot(projection=aia_sequence_aligned[0])
ani = aia_sequence_aligned.plot(axes=ax, cmap='sdoaia171', norm=ImageNormalize(vmin=0, vmax=5e3, stretch=SqrtStretch()))

plt.show()
