"""
===========================
Adding an Earth scale image
===========================

This example shows how to plot a map with an image of the Earth added for scale.
"""
from urllib.request import urlretrieve

import matplotlib.pyplot as plt
from matplotlib.image import AxesImage
from PIL import Image

import astropy.units as u
from astropy.constants import R_earth
from astropy.coordinates import SkyCoord
from astropy.wcs import WCS

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
# Now we plot the AIA image and superimpose the Earth for scale. We manually
# create and add the image artist so that the Earth will not be considered
# when autoscaling plot limits. We also add a text label above the Earth.

fig = plt.figure()
ax = fig.add_subplot(projection=cutout_map)
cutout_map.plot(clip_interval=(1, 99.9)*u.percent)

earth_artist = AxesImage(ax, transform=ax.get_transform(WCS(earth_wcs)))
earth_artist.set_data(earth)
ax.add_artist(earth_artist)

ax.text(
    earth_x.to_value('deg'), (earth_y + earth_diameter).to_value('deg'),
    'Earth to scale', color='white', fontsize=12, horizontalalignment='center',
    transform=ax.get_transform('world')
)

plt.show()
