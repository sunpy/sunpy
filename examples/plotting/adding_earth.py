"""
===========================
Adding an Earth scale image
===========================

This example shows how to plot a map with an image of the Earth added for scale.

This example was adapted from an example written by Alex Russell.
"""
from urllib.request import urlretrieve

import matplotlib.pyplot as plt
import numpy as np
from PIL import Image

import astropy.units as u
from astropy.constants import R_earth
from astropy.coordinates import SkyCoord

import sunpy.data.sample
import sunpy.map
from sunpy.coordinates.utils import solar_angle_equivalency

###############################################################################
# We start with a sample AIA image.

cutout_map = sunpy.map.Map(sunpy.data.sample.AIA_193_CUTOUT01_IMAGE)

###############################################################################
# We download a low-resolution Earth image from Wikimedia Commons that already
# has a transparent background. We also crop the image tightly.

earth_file = urlretrieve("https://upload.wikimedia.org/wikipedia/commons/thumb/4/43/The_Earth_seen_from_Apollo_17_with_transparent_background.png/250px-The_Earth_seen_from_Apollo_17_with_transparent_background.png")[0]
earth = Image.open(earth_file)
earth = earth.crop(earth.getbbox())

##############################################################################
# The first step in plotting the Earth is to convert the Earth's diameter in km
# to arcsec on the plane-of-sky for the time of the observation.  We use the
# observer coordinate from the map. For more information about the
# transformation, see :func:`~sunpy.coordinates.utils.solar_angle_equivalency`.

earth_diameter = 2 * R_earth.to_value(u.arcsec, equivalencies=solar_angle_equivalency(cutout_map.observer_coordinate))
print(f"Earth diameter = {earth_diameter:.2f} arcsec")

##############################################################################
# Now we plot the AIA image and superimpose the Earth for scale. The image of
# the Earth is plotted by inserting a new axis within the AIA figure.
#
# The width and height of the new axis is specified in WCS data coordinates,
# which must be calculated to agree with the Earth's diameter.
# This example does the calculation by the following steps.
#
# 1. An invisible quadrangle is drawn with the correct size and position.
#    Its extents provide the required display size of Earth.
#
# 2. The image of the Earth is plotted on a new insert axis within the AIA map.
#    `The tricky part is navigating the different coordinate systems used by matplotlib
#    and transforming between them. <https://matplotlib.org/stable/users/explain/artists/transforms_tutorial.html>`__

# Specify arcsec location to display the Earth image
earth_x = 1000
earth_y = -200

# Plot AIA
fig = plt.figure()
ax = plt.axes(projection=cutout_map)
cutout_map.plot(clip_interval=(1, 99.9)*u.percent)

# Create new axis of correct size to display Earth image
coords = SkyCoord(Tx=(earth_x, earth_x + earth_diameter) * u.arcsec, Ty=(earth_y, earth_y + earth_diameter) * u.arcsec, frame=cutout_map.coordinate_frame)
rect = cutout_map.draw_quadrangle(coords, axes=ax, lw=0)
fc = rect.get_extents().transformed(ax.transData.inverted()).corners()
earth_ax = ax.inset_axes(fc[0].tolist() + (fc[-1]-fc[0]).tolist(), transform=ax.transData)

# Now we want to add the image of the Earth
earth_ax.imshow(earth)
earth_ax.axis('off')
earth_ax.set_title('Earth to scale', color='white', fontsize=12)

plt.show()
