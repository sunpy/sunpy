"""
===========================
Adding an Earth scale image
===========================

This example shows how to plot a map with an image of the Earth added for scale.

The Earth image used in this example is produced by NASA's EPIC Team.
Their images are freely available for re-production or re-use, including commercial purposes.
The NASA EPIC Team ask that they be be given credit for the original materials.
`For information about usage. <https://epic.gsfc.nasa.gov/about>`__

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
from sunpy.coordinates import get_body_heliographic_stonyhurst
from sunpy.coordinates.utils import solar_angle_equivalency

###############################################################################
# We start with the sample data and downloading the Earth image.

cutout_map = sunpy.map.Map(sunpy.data.sample.AIA_193_CUTOUT01_IMAGE)
# Here we just pick a nice recent(ish) image.
urlretrieve(
    url="https://epic.gsfc.nasa.gov/archive/enhanced/2025/06/17/png/epic_RGB_20250617102539.png",
    filename="epic_RGB_20250617102539.png"
)

##############################################################################
# The first step in plotting the Earth is to convert the Earth's diameter in km
# to arcsec on the plane-of-sky for the time of the observation.
#
# This is ~17.65 arcsec for this example, but varies according to the Sun-Earth
# distance, so it needs to be evaluated for the observation date. For more information
# about the transformation, see `~sunpy.coordinates.utils.solar_angle_equivalency`.

observer = get_body_heliographic_stonyhurst('earth', cutout_map.meta['date'])
earth_diameter = (R_earth * 2).to(u.arcsec, equivalencies=solar_angle_equivalency(observer)).value

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
# 2. The image of the Earth is plotted on the new insert axis on top of the AIA map.
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

# Plot the image of the Earth
img = np.asarray(Image.open("epic_RGB_20250617102539.png"))
# Add an alpha channel and set that to be 0 where there is no data
img = np.concatenate([img, np.where(np.sum(img, axis=-1) == 0, 0, 255)[..., None]], axis=-1)
x_extent  = np.argwhere(np.sum(img[:,img.shape[1]//2, ...], axis=1)).flatten()[[0, -1]]
y_extent  = np.argwhere(np.sum(img[shape[0]//2, ...], axis=1)).flatten()[[0, -1]]
img_cut = img[y_extent[0]:y_extent[1], x_extent[0]:x_extent[1]]
earth_ax.imshow(img_cut)
earth_ax.axis('off')

# Here the location was hand picked to be in the middle of the Earth image.
earth_ax.annotate('Earth to scale', [-2.1, 1.2], xycoords = 'axes fraction', color='white', fontsize=12, transform=earth_ax.transAxes)

plt.show()
