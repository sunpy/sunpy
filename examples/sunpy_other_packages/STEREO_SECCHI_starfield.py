# -*- coding: utf-8 -*-
"""
===========================================================
Identifying stars in a STEREO/SECCHI COR2 coronagraph image
===========================================================

How to use the `Vizier star catalog <https://vizier.u-strasbg.fr>`_ to
identify bright points in a STEREO/SECCHI COR 2 image as stars.
As a bonus, we also identify Mars.
"""
import matplotlib.pyplot as plt

# astroquery is not a dependency of SunPy so will need to be
# install this package separately.
from astroquery.vizier import Vizier
import astropy.units as u
from astropy.coordinates import SkyCoord, get_body_barycentric

import sunpy.map
from sunpy.coordinates import get_body_heliographic_stonyhurst, frames
from sunpy.net import helioviewer

###############################################################################
# Let's download a STEREO-A SECCHI COR2 image from Helioviewer which provide
# pre-processed images and load it into a Map.
hv = helioviewer.HelioviewerClient()
f = hv.download_jp2('2014/05/15 07:54', observatory='STEREO_A', instrument='SECCHI', detector='COR2')
cor2 = sunpy.map.Map(f)

###############################################################################
# To efficiently search the star field we need to know what stars are near the
# Sun as observed by STEREO. We need the vector that points from STEREO to the Sun.
# By converting to HCRS we get the vector from the Sun to STEREO
sun_to_stereo = cor2.observer_coordinate.transform_to('hcrs')

###############################################################################
# We next reflect the vector to get our search vector which points from STEREO
# to the Sun
stereo_to_sun = SkyCoord(-sun_to_stereo.data, obstime=sun_to_stereo.obstime, frame='hcrs')

###############################################################################
# Let's look up bright stars using the Vizier search capability provided by
# astroquery (note that this is not a required package of SunPy so you will likely
# need to install it). We will search the GAIA2 star catalog for stars with magnitude
# brighter than 7.
vv = Vizier(columns=['**'], row_limit=-1, column_filters={'Gmag': '<7'}, timeout=1200)
vv.ROW_LIMIT = -1
result = vv.query_region(stereo_to_sun, radius=4 * u.deg, catalog='I/345/gaia2')

###############################################################################
# Let's see how many stars we've found.
print(len(result[0]))

###############################################################################
# Now we load each star into a coordinate and transform it into the COR2
# image coordinates. Since we don't know the distance to each of these stars
# we will just put them very far away.
hpc_coords = []
for this_object in result[0]:
    tbl_crds = SkyCoord(this_object['RA_ICRS'] * u.deg, this_object['DE_ICRS'] * u.deg,
                        1e12 * u.km, frame='icrs', obstime=cor2.date)
    hpc_coords.append(tbl_crds.transform_to(cor2.coordinate_frame))

###############################################################################
# One of the bright features is actually Mars so let's also get that coordinate.
# get the location of Mars.
mars = get_body_heliographic_stonyhurst('mars', cor2.date, observer=cor2.observer_coordinate)
mars_hpc = mars.transform_to(frames.Helioprojective(observer=cor2.observer_coordinate))

###############################################################################
# Let's plot the results.
ax = plt.subplot(projection=cor2)

# Let's tweak the axis to show in degrees instead of arcsec
lon, lat = ax.coords
lon.set_major_formatter('d.dd')
lat.set_major_formatter('d.dd')

cor2.plot(axes=ax, vmin=0, vmax=600)
cor2.draw_limb()

# plot the position of Mars
ax.plot_coord(mars_hpc, 's', color='white', fillstyle='none', markersize=12, label='Mars')
# Plot all of the stars
for this_coord in hpc_coords:
    ax.plot_coord(this_coord, 'o', color='white', fillstyle='none')
plt.legend()
plt.show()
