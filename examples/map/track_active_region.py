"""
===============================================
Tracking an Active Region Across the Solar Disk
===============================================

This example demonstrates how to track an active region as it rotates across the solar disk
and make cutouts around that active region at each time step to build a tracked datacube.
"""
# sphinx_gallery_thumbnail_number = 3
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
# of half of a day at a cadence of 1 image every 1 hour.

query = Fido.search(a.Time('2018-05-30 00:00:00', '2018-05-30 12:00:00'),
                    a.Instrument.aia,
                    a.Wavelength(171*u.angstrom),
                    a.Sample(1*u.h))
print(query)
files = Fido.fetch(query)

###############################################################################
# Now that we have a set of images in time, we can create a `~sunpy.map.MapSequence` to hold all of them
# and animate that sequence in time.

aia_sequence = sunpy.map.Map(files, sequence=True)

fig = plt.figure()
ax = fig.add_subplot(projection=aia_sequence[0])
norm = norm=ImageNormalize(vmin=0, vmax=3e3, stretch=SqrtStretch())
ani = aia_sequence.plot(axes=ax, norm=norm)

###############################################################################
# Next, let's crop one of the maps in our sequence to the active region of interest.

corner = SkyCoord(Tx=-250*u.arcsec, Ty=0*u.arcsec, frame=aia_sequence[6].coordinate_frame)
cutout_map = aia_sequence[6].submap(corner, width=500*u.arcsec, height=500*u.arcsec)

fig = plt.figure()
ax = fig.add_subplot(projection=cutout_map)
cutout_map.plot(axes=ax)

###############################################################################
# Now that we have our cutout around our active region, we can reproject each map in our sequence to
# the WCS of that cutout. Additionally, we will use the `~sunpy.coordinates.propagate_with_solar_surface`
# context manager to adjust the field of view of the cutout with the rotation of the solar surface.
# We use the ``preserve_date_obs`` keyword argument to preserve the original observation time in each
# reprojected map. Otherwise, the observation time of the reprojected maps would all be the same as the
# cutout WCS. Note that using this keyword does not affect the resulting coordinate frame of the reprojected
# map which will be defined by the date of the cutout WCS.

with propagate_with_solar_surface():
    aia_sequence_aligned = sunpy.map.Map(
        [m.reproject_to(cutout_map.wcs, preserve_date_obs=True) for m in aia_sequence],
        sequence=True
    )

###############################################################################
# Finally, we can animate our sequence of reprojected cutouts to confirm that we've tracked our active
# region of interest as a function of time across the solar disk.

fig = plt.figure()
ax = fig.add_subplot(projection=aia_sequence_aligned[0])
ani = aia_sequence_aligned.plot(axes=ax, cmap='sdoaia171', norm=norm)

plt.show()
