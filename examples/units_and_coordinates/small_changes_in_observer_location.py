"""
============================================
Handling Small Changes in Observer Locations
============================================

In this example we will investigate the effects of small changes in observer locations when performing coordinate transforms.
We will demonstrate these effects and then describe some different assumptions that can be made to make working with data from different observers easier.
"""
# sphinx_gallery_tags = ["Coordinates", "AIA", "Map"]

import astropy.units as u
from astropy.coordinates import SkyCoord

import sunpy.map
from sunpy.data.sample import (
    AIA_094_IMAGE,
    AIA_131_IMAGE,
    AIA_171_IMAGE,
    AIA_193_IMAGE,
    AIA_211_IMAGE,
    AIA_304_IMAGE,
    AIA_335_IMAGE,
)

################################################################################
# For our example data, we will use a collection of AIA passbands
# taken close together in time.

aia_channels = [
    AIA_094_IMAGE,
    AIA_131_IMAGE,
    AIA_171_IMAGE,
    AIA_193_IMAGE,
    AIA_211_IMAGE,
    AIA_304_IMAGE,
    AIA_335_IMAGE,
]
aia_channel_maps = [sunpy.map.Map(f) for f in aia_channels]
aia_channel_maps = list(sorted(aia_channel_maps, key=lambda smap: smap.date))

for smap in aia_channel_maps:
    print(f"AIA {smap.wavelength.value:3.0f} {smap.wavelength.unit} has Time: {smap.date}")

time_delta = (aia_channel_maps[-1].date - aia_channel_maps[0].date).to(u.s)
print(f"Time difference between first and last map: {time_delta}")

################################################################################
# Now, let's extract the first and last observer coordinates and consider how much they change.

first_observer = aia_channel_maps[0].observer_coordinate
last_observer = aia_channel_maps[-1].observer_coordinate
separation = first_observer.separation_3d(last_observer)

print(f"Physical separation between first and last map: {separation.to(u.km)}")
print(f"SDO's speed {(separation / time_delta).to(u.km / u.s)}")

################################################################################
# Just over 1 Mm seems like a large separation, but the angular
# separation as seen by the center of the Sun is still very small:

angular_separation = first_observer.separation(last_observer).to(u.arcsec)
print(f"Angular separation between the two observers (from solar center): {angular_separation}")


hgs_0 = SkyCoord(
    0 * u.deg,
    0 * u.deg,
    frame="heliographic_stonyhurst",
    obstime=aia_channel_maps[0].reference_date,
)

hpc_first = hgs_0.transform_to(aia_channel_maps[0].coordinate_frame)
hpc_last = hgs_0.transform_to(aia_channel_maps[-1].coordinate_frame)
