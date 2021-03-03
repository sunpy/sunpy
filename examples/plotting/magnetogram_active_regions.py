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

if 'I' in srs_table['ID'] or 'IA' in srs_table['ID']:
    srs_table = srs_table[np.logical_or(srs_table['ID'] == 'I',
                                        srs_table['ID'] == 'IA')]
else:
    raise ValueError("No I or IA entries for this date.")

##############################################################################
# Now we extract the latitudes, longitudes and the region numbers. We make an
# empty list if there are no ARs.

if srs_table is not None:
    lats = srs_table['Latitude']
    lngs = srs_table['Longitude']
    numbers = srs_table['Number']
else:
    lats = lngs = numbers = []

##############################################################################
# Let's plot the results by defining coordinates for each location.

ax = plt.subplot(projection=smap)
smap.plot(vmin=-120, vmax=120)
smap.draw_limb()
ax.set_autoscale_on(False)

if len(lats) > 0:
    c = SkyCoord(lngs, lats, frame="heliographic_stonyhurst")
    ax.plot_coord(c, 'o')

    for i, num in enumerate(numbers):
        ax.annotate(num, (lngs[i].value, lats[i].value),
                    xycoords=ax.get_transform('heliographic_stonyhurst'),
                    color='red',
                    fontweight='bold',
                    )
plt.show()
