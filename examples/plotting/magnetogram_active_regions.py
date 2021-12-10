"""
==========================================================
Overplotting SRS active region locations on a magnetograms
==========================================================

How to find and plot the location of an active region on an HMI magnetogram.
"""
import matplotlib.pyplot as plt
import numpy as np

from astropy.coordinates import SkyCoord

import sunpy.coordinates
import sunpy.data.sample
import sunpy.map
from sunpy.io.special import srs

##############################################################################
# For this example, we will start with the sample data. We need an HMI file and
# use it to create a map, and the SRS table which contains a list of active
# regions. Both of these data can be downloaded with ``Fido``.

smap = sunpy.map.Map(sunpy.data.sample.HMI_LOS_IMAGE)
srs_table = srs.read_srs(sunpy.data.sample.SRS_TABLE)

##############################################################################
# We only need the rows which have 'ID' = 'I' or 'IA'.
# Some tables do not have these columns, so we exit the script if they are not
# present.

srs_table = srs_table[np.logical_or(srs_table['ID'] == 'I', srs_table['ID'] == 'IA')]

##############################################################################
# Now we extract the latitudes, longitudes and the region numbers.
# :meth:`~astropy.visualization.wcsaxes.WCSAxes.plot_coord` will error on some
# versions of Astropy if a coordinate component contains a mask, but since
# none of the masks in these arrays here actually mask any elements, we simply
# remove the masks.

lats = srs_table['Latitude']
if hasattr(lats, 'mask'):
    lats = lats.unmasked
lngs = srs_table['Longitude']
if hasattr(lngs, 'mask'):
    lngs = lngs.unmasked
numbers = srs_table['Number']

##############################################################################
# Let's plot the results by defining coordinates for each location.

ax = plt.subplot(projection=smap)
smap.plot(vmin=-120, vmax=120)
smap.draw_limb()
c = SkyCoord(lngs, lats, frame="heliographic_stonyhurst")
ax.plot_coord(c, 'o')
for num, lng, lat in zip(numbers, lngs.value, lats.value):
    ax.annotate(num, (lng, lat),
                xycoords=ax.get_transform('heliographic_stonyhurst'),
                color='red',
                fontweight='bold',
                )

plt.show()
