"""
===========================
Adding an Earth scale image
===========================

This example shows how to plot a map with an image of the Earth added for scale.
"""
import matplotlib.pyplot as plt
import numpy as np
from PIL import Image

import astropy.units as u
from astropy.constants import R_earth
from astropy.coordinates import SkyCoord
from astropy.wcs import WCS

import sunpy.map
from sunpy.coordinates.utils import solar_angle_equivalency
from sunpy.data import EARTH_IMAGE
from sunpy.data.sample import AIA_193_CUTOUT01_IMAGE

###############################################################################
# We start with a sample AIA image.

cutout_map = sunpy.map.Map(AIA_193_CUTOUT01_IMAGE)

###############################################################################
# We use a (low-resolution) image of Earth that we provide with `sunpy`. You
# can use a different image as desired, and it should have a transparent
# background. We also crop the image tightly so that we can assume that the
# image width/height are equal to the Earth diameter. Finally, we flip the
# image vertically so that it is oriented correctly when plotted with the pixel
# origin in the lower-left corner, which is the convention for maps.

earth = Image.open(EARTH_IMAGE)
earth = earth.crop(earth.getbbox())
earth = earth.transpose(Image.FLIP_TOP_BOTTOM)

##############################################################################
# The first step in plotting the Earth is to convert the Earth's diameter in km
# to arcsec on the plane-of-sky for the time of the observation.  We use the
# observer coordinate from the map. For more information about the
# transformation, see :func:`~sunpy.coordinates.utils.solar_angle_equivalency`.

earth_diameter = 2 * R_earth.to(u.arcsec, equivalencies=solar_angle_equivalency(cutout_map.observer_coordinate))
print(f"Earth diameter = {earth_diameter:.2f}")

##############################################################################
# We then create a WCS header so that the Earth will be plotted correctly
# on any underlying map projection. We center the Earth on a specific point in
# the map.

earth_x = 1000 * u.arcsec
earth_y = -200 * u.arcsec
earth_center = SkyCoord(earth_x, earth_y, frame=cutout_map.coordinate_frame)
earth_wcs = sunpy.map.make_fitswcs_header(
    (earth.height, earth.width), earth_center,
    scale=earth_diameter / ([earth.width, earth.height] * u.pix)
)

##############################################################################
# Now we plot the AIA image and superimpose the Earth for scale. We use
# :meth:`matplotlib.axes.Axes.pcolormesh` to plot the pixels of the Earth
# image as a quadrilateral mesh, which accommodates any WCS of the underlying
# map. We also add a text label above the Earth.

fig = plt.figure()
ax = fig.add_subplot(projection=cutout_map)
cutout_map.plot(clip_interval=(1, 99.9)*u.percent)

ax.pcolormesh(
    *np.meshgrid(np.arange(earth.width + 1) - 0.5, np.arange(earth.height + 1) - 0.5),
    earth, shading='flat', transform=ax.get_transform(WCS(earth_wcs))
)

ax.text(
    earth_x.to_value('deg'), (earth_y + earth_diameter).to_value('deg'),
    'Earth to scale', color='white', fontsize=12, horizontalalignment='center',
    transform=ax.get_transform('world')
)

plt.show()
