"""
====================================
Making a cutout with the HMI instrument
====================================

How to use the `~sunpy.net.jsoc.JSOCClient` to make a cutout from HMI data.
"""
# sphinx_gallery_tags = ["Cutout", "HMI", "JSOC", "Download"]

import matplotlib.pyplot as plt

import astropy.units as u
from astropy.coordinates import SkyCoord

import sunpy.map
import sunpy.physics.differential_rotation as drot
from sunpy.data.sample import HMI_MAGNETOGRAM_IMAGE
from sunpy.net import Fido
from sunpy.net import attrs as a

###############################################################################
# We can use the sample data to plot a cutout region.
# Let's first create a map of the sample data.

magnetogram = sunpy.map.Map(HMI_MAGNETOGRAM_IMAGE)

###############################################################################
# Now let's create a cutout.  We can do this with the
# `~sunpy.net.jsoc.JSOCClient` class.

right_corner = SkyCoord(Tx=158*u.arcsec, Ty=350*u.arcsec, frame=magnetogram.coordinate_frame)

###############################################################################
# Let's create the cutout request and download the data.

hpc_coords = sunpy.map.all_coordinates_from_map(magnetogram)
mask = ~sunpy.coordinates.utils.coordinate_is_on_solar_disk(hpc_coords)
magnetogram_big = sunpy.map.Map(magnetogram.data, magnetogram.meta, mask=mask)

###############################################################################
# Let's plot the results.

fig = plt.figure()
ax = fig.add_subplot(projection=magnetogram_big)
magnetogram_big.plot(axes=ax, title="HMI Magnetogram with Mask")
plt.show()
