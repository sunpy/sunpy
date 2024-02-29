"""
=======================================================
Setting the correct position for SOHO in a LASCO C3 Map
=======================================================

How to get the correct location of SOHO using JPL HORIZONS
and update the header.
"""
# sphinx_gallery_thumbnail_number = 2

import hvpy
import matplotlib.pyplot as plt
import numpy as np

import sunpy.map
from sunpy.coordinates.ephemeris import get_body_heliographic_stonyhurst, get_horizons_coord
from sunpy.time import parse_time
from sunpy.util.config import get_and_create_download_dir

###############################################################################
# Let's download a SOHO/LASCO C3 image Helioviewer.org and load it into a Map.
# The reason to use Helioviewer.org is that they provide processed images.
# We download to the default sunpy download directory.

lasco_file = hvpy.save_file(hvpy.getJP2Image(parse_time('2000/02/27 07:42').datetime,
                                             hvpy.DataSource.LASCO_C3.value),
                            get_and_create_download_dir() + "/LASCO_C3.jp2")
lasco_map = sunpy.map.Map(lasco_file)

###############################################################################
# A user warning let's you know that there is missing metadata for the observer
# location. sunpy goes ahead and assumes that the observer is at Earth.

print(lasco_map.observer_coordinate)

###############################################################################
# To check that this worked let's get the location of Mercury in this exposure
# and show that it is in the correct location.

mercury_wrong = get_body_heliographic_stonyhurst(
    'mercury', lasco_map.date, observer=lasco_map.observer_coordinate)
mercury_hpc_wrong = mercury_wrong.transform_to(lasco_map.coordinate_frame)
print(mercury_hpc_wrong)

##############################################################################
# Let's plot how this looks with the incorrect observer information.

fig = plt.figure()
ax = fig.add_subplot(projection=lasco_map)

# Let's tweak the axis to show in degrees instead of arcsec
lon, lat = ax.coords
lon.set_major_formatter('d.dd')
lat.set_major_formatter('d.dd')
ax.plot_coord(mercury_hpc_wrong, 's', color='white',
              fillstyle='none', markersize=12, label='Mercury')
lasco_map.plot(axes=ax)

plt.show()

###############################################################################
# SOHO is actually a `halo orbit <https://en.wikipedia.org/wiki/Solar_and_Heliospheric_Observatory#Orbit>`__
# around the Sun–Earth L1 point, about 1 million km away from the Earth.
# The following functions queries JPL HORIZONS which includes positions of major spacecraft.
# This function requires an internet connection to fetch the ephemeris data.

soho = get_horizons_coord('SOHO', lasco_map.date)

###############################################################################
# For fun, let's see how far away from the Earth SOHO is by converting to
# an Earth-centered coordinate system (GCRS).

print(soho.transform_to('gcrs').distance.to('km'))

###############################################################################
# Let's fix the header.

lasco_map.meta['HGLN_OBS'] = soho.lon.to('deg').value
lasco_map.meta['HGLT_OBS'] = soho.lat.to('deg').value
lasco_map.meta['DSUN_OBS'] = soho.radius.to('m').value

###############################################################################
# Let's get the right position now.

mercury = get_body_heliographic_stonyhurst(
    'mercury', lasco_map.date, observer=lasco_map.observer_coordinate)
mercury_hpc = mercury.transform_to(lasco_map.coordinate_frame)

###############################################################################
# The difference between the incorrect position and the right one is:

r = np.sqrt((mercury_hpc.Tx - mercury_hpc_wrong.Tx) ** 2 +
            (mercury_hpc.Ty - mercury_hpc_wrong.Ty) ** 2)
print(r)

##############################################################################
# Let's plot the results.

fig = plt.figure()
ax = fig.add_subplot(projection=lasco_map)

# Let's tweak the axis to show in degrees instead of arcsec
lon, lat = ax.coords
lon.set_major_formatter('d.dd')
lat.set_major_formatter('d.dd')
ax.plot_coord(mercury_hpc, 's', color='white', fillstyle='none', markersize=12, label='Mercury')
lasco_map.plot(axes=ax)

plt.show()
