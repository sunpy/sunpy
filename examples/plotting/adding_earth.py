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
# We start with the sample data and download a low-resolution Earth image from
# Wikimedia Commons.

cutout_map = sunpy.map.Map(sunpy.data.sample.AIA_193_CUTOUT01_IMAGE)
earth_image = urlretrieve("https://upload.wikimedia.org/wikipedia/commons/thumb/7/7e/Earth_Western_Hemisphere_2002.png/250px-Earth_Western_Hemisphere_2002.png")[0]

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
earth_img = np.asarray(Image.open(earth_image))
# Add an alpha channel and set that to be 0 where there is no data,
# this removes the background
earth_img = np.concatenate([earth_img, np.where(np.sum(earth_img, axis=-1) == 0, 0, 255)[..., None]], axis=-1)
# Now we crop the image just to save on computation
x_extent  = np.argwhere(np.sum(earth_img[:,earth_img.shape[1]//2, ...], axis=1)).flatten()[[0, -1]]
y_extent  = np.argwhere(np.sum(earth_img[earth_img.shape[0]//2, ...], axis=1)).flatten()[[0, -1]]
earth_crop = earth_img[y_extent[0]:y_extent[1], x_extent[0]:x_extent[1]]

earth_ax.imshow(earth_crop)
earth_ax.axis('off')
earth_ax.set_title('Earth to scale', color='white', fontsize=12)

plt.show()
