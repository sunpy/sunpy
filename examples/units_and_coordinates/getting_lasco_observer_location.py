# -*- coding: utf-8 -*-
"""
=======================================================
Setting the correct position for SOHO in a LASCO C3 Map
=======================================================

How to get the correct location of SOHO using JPL HORIZONS
and update the header.

This requires the intallation of the
`astroquery <https://astroquery.readthedocs.io/en/latest/>`__
package and an internet connection.
`astroquery <https://astroquery.readthedocs.io/en/latest/>`__ can be installed ontop of
the existing sunpy conda environemnt: conda install -c astropy astroquery
"""

# sphinx_gallery_thumbnail_number = 2
import numpy as np
import matplotlib.pyplot as plt

import sunpy.map
from sunpy.coordinates.ephemeris import get_horizons_coord, get_body_heliographic_stonyhurst
from sunpy.net import helioviewer
hv = helioviewer.HelioviewerClient()

###############################################################################
# Let's download a SOHO/LASCO C3 image from Helioviewer which provided
# pre-processed images and load it into a Map.
f = hv.download_jp2('2000/02/27 07:42', observatory='SOHO', instrument='LASCO', detector='C3')
lasco = sunpy.map.Map(f)

###############################################################################
# A user warning let's you know that there is missing metadata for the observer
# location. SunPy goes ahead and assumes that the observer is at Earth.
print(lasco.observer_coordinate)

###############################################################################
# To check that this worked let's get the location of Mercury in this exposure
# and show that it is in the correct location.
mercury_wrong = get_body_heliographic_stonyhurst('mercury', lasco.date, observer=lasco.observer_coordinate)
mercury_hpc_wrong = mercury_wrong.transform_to(lasco.coordinate_frame)
print(mercury_hpc_wrong)

##############################################################################
# Let's plot how this looks with the incorrect observer information.
ax = plt.subplot(projection=lasco)

# Let's tweak the axis to show in degrees instead of arcsec
lon, lat = ax.coords
lon.set_major_formatter('d.dd')
lat.set_major_formatter('d.dd')
ax.plot_coord(mercury_hpc_wrong, 's', color='white', fillstyle='none', markersize=12, label='Mercury')
lasco.plot()
plt.show()

###############################################################################
# SOHO is actually a `halo orbit <https://en.wikipedia.org/wiki/Solar_and_Heliospheric_Observatory#Orbit>`__
# around the Sun–Earth L1 point, about 1 million km away from the Earth.
# The following functions queries JPL HORIZONS which includes positions of major spacecraft.
# This function requires an internet connection to fetch the ephemeris data.
soho = get_horizons_coord('SOHO', lasco.date)

###############################################################################
# For fun, let's see how far away from the Earth SOHO is by converting to
# an Earth-centered coordinate system (GCRS).
print(soho.transform_to('gcrs').distance.to('km'))

###############################################################################
# Let's fix the header.
lasco.meta['HGLN_OBS'] = soho.lon.to('deg').value
lasco.meta['HGLT_OBS'] = soho.lat.to('deg').value
lasco.meta['DSUN_OBS'] = soho.radius.to('m').value

###############################################################################
# Let's get the right position now.
mercury = get_body_heliographic_stonyhurst('mercury', lasco.date, observer=lasco.observer_coordinate)
mercury_hpc = mercury.transform_to(lasco.coordinate_frame)

###############################################################################
# The difference between the incorrect position and the right one is:
r = np.sqrt((mercury_hpc.Tx - mercury_hpc_wrong.Tx) ** 2 + (mercury_hpc.Ty - mercury_hpc_wrong.Ty) ** 2)
print(r)

##############################################################################
# Let's plot the results.
ax = plt.subplot(projection=lasco)

# Let's tweak the axis to show in degrees instead of arcsec
lon, lat = ax.coords
lon.set_major_formatter('d.dd')
lat.set_major_formatter('d.dd')
ax.plot_coord(mercury_hpc, 's', color='white', fillstyle='none', markersize=12, label='Mercury')
lasco.plot()
plt.show()
