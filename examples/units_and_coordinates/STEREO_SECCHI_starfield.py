"""
===========================================================
Identifying stars in a STEREO/SECCHI COR2 coronagraph image
===========================================================

Since the field of view in SECCHI COR2 images can span 2 to 15 solar radii,
we often observe stars in these image data. In this example, we will use the
`Astroquery package <https://astroquery.readthedocs.io/en/latest/>`__ to query the
`VizieR star catalog <http://vizier.u-strasbg.fr/viz-bin/VizieR>`__ for stars observed
by the `Gaia satellite <https://sci.esa.int/web/gaia/>`__ within the SECCHI COR2 field of view.
Then we will use the coordinates framework in SunPy and AstroPy to transform the coordinates
returned by VizieR into SECCHI COR2 image coordinates. As a bonus, we'll also identify Mars."

This requires the installation of the `astroquery <https://astroquery.readthedocs.io/en/latest/>`__
package, which can be installed on top of the existing sunpy conda
environment: ``conda install -c astropy astroquery`` and an active internet connection.
"""
import matplotlib.pyplot as plt
from astroquery.vizier import Vizier

import astropy.units as u
from astropy.coordinates import SkyCoord

import sunpy.map
from sunpy.coordinates import frames, get_body_heliographic_stonyhurst
from sunpy.net import helioviewer

###############################################################################
# Let's download a STEREO-A SECCHI COR2 image from Helioviewer which provide
# pre-processed images and load it into a Map.

hv = helioviewer.HelioviewerClient()
f = hv.download_jp2('2014/05/15 07:54', observatory='STEREO_A',
                    instrument='SECCHI', detector='COR2')
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

# Plot the position of Mars
ax.plot_coord(mars_hpc, 's', color='white', fillstyle='none', markersize=12, label='Mars')
# Plot all of the stars
for this_coord in hpc_coords:
    ax.plot_coord(this_coord, 'o', color='white', fillstyle='none')
plt.legend()
plt.show()
