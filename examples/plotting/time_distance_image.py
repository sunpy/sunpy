"""
=====================================================
Creating a time-distance plot from a sequence of maps
=====================================================

This example showcases how you can use :func:`sunpy.map.pixelate_coord_path` and
:func:`sunpy.map.sample_at_coords()` on a sequence of images to create a time-distance diagram.
"""
import matplotlib.pyplot as plt
import numpy as np

import astropy.units as u
from astropy.coordinates import SkyCoord

import sunpy.map
from sunpy.coordinates import propagate_with_solar_surface
from sunpy.coordinates.screens import SphericalScreen
from sunpy.net import Fido
from sunpy.net import attrs as a

###############################################################################
# First, we will need to acquire a sequence of images from SDO/AIA.

query = Fido.search(
    a.Time('2025-05-01 00:00:00', '2025-05-01 00:10:00'),
    a.Instrument.aia,
    a.Wavelength(171*u.angstrom)
)
files = Fido.fetch(query)
aia_sequence = sunpy.map.Map(files, sequence=True)

###############################################################################
# Our target will be a set of loops in the corona over a small bright region.
# We will create a rectangular slice (i.e., a map cutout) around this region.

corner = SkyCoord(Tx=800*u.arcsec, Ty=200*u.arcsec, frame=aia_sequence[0].coordinate_frame)
cut_sequence = []
for aia_map in aia_sequence:
    with SphericalScreen(aia_map.observer_coordinate):
        cut_sequence.append(aia_map.submap(corner, width=400*u.arcsec, height=400*u.arcsec))
cut_aia_sequence = sunpy.map.Map(cut_sequence, sequence=True)

###############################################################################
# Now that we have our cutouts, we can reproject each map in our sequence to
# the WCS of that cutout using the `~sunpy.coordinates.propagate_with_solar_surface`
# context manager to adjust the field of view of the cutout with the rotation of the solar surface.

with propagate_with_solar_surface():
    with SphericalScreen(cut_aia_sequence[0].observer_coordinate):
        aia_sequence_aligned = sunpy.map.Map([m.reproject_to(cut_aia_sequence[0].wcs) for m in cut_aia_sequence], sequence=True)

###############################################################################
# Next, we will define a path in space and then use that to fetch the pixel
# coordinates for every pixel that intersects with the physical path.

# Note this is [X0, X1] and [Y0, Y1] in the coordinate frame of the map.
line_coords = SkyCoord([900, 1000],[456, 250] , unit=u.arcsec, frame=aia_sequence_aligned[0].coordinate_frame)

###############################################################################
# We can now plot ``line_coords`` on the first map to see what we are intersecting.

fig = plt.figure(figsize=(10, 4))
ax1 = fig.add_subplot(111, projection=aia_sequence_aligned[0])
aia_sequence_aligned[0].plot(axes=ax1)
ax1.plot_coord(line_coords)

###############################################################################
# Now let us construct the time-distance diagram. First we need to get the data
# over the path we defined above for the entire sequence.

intensities = []
for aia_map in aia_sequence_aligned:
    intensity_coords = sunpy.map.pixelate_coord_path(aia_map, line_coords)
    intensities.append(sunpy.map.sample_at_coords(aia_map, intensity_coords))

# The above will give us a list of 1D arrays, one for each map in the sequence.
# We need to stack them into a 2D array.
intensities = np.vstack(intensities)

# This defines the distance along the path in pixels.
angular_separation = intensity_coords.separation(intensity_coords[0]).to(u.arcsec)

###############################################################################
# We can now plot the time-distance diagram.

fig = plt.figure(figsize=(10, 4))
ax1 = fig.add_subplot(111, )
ax1.imshow(intensities, aspect='auto', origin='lower', extent=[0, len(aia_sequence_aligned) * 12*u.s, 0, angular_separation[-1].value])
ax1.set_xlabel('Time (s)')
ax1.set_ylabel('Distance (arcsec)')

plt.show()
