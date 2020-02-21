"""
==========================================================
Overplotting SRS active region locations on a magnetograms
==========================================================

How to find and plot the location of an active region on an HMI magnetogram.
"""
import numpy as np
import matplotlib.pyplot as plt

import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.time import TimeDelta

import sunpy.map
import sunpy.coordinates
from sunpy.io.special import srs
from sunpy.time import parse_time
from sunpy.net import Fido, attrs as a

##############################################################################
# For this example, we will search for and download a single HMI using Fido.
start_time = parse_time("2017-01-25")
end_time = start_time + TimeDelta(23*u.hour + 59*u.minute + 59*u.second)
results = Fido.search(a.Time(start_time, end_time),
                      a.Instrument('HMI') & a.vso.Physobs("LOS_magnetic_field"),
                      a.Sample(60 * u.second))

##############################################################################
# Let's select only the first file, download it and create a map.
result = results[0, 0]
file_name = Fido.fetch(result)
smap = sunpy.map.Map(file_name)

##############################################################################
# Download the SRS file.
srs_results = Fido.search(a.Time(start_time, end_time), a.Instrument('SRS_TABLE'))
srs_downloaded_files = Fido.fetch(srs_results)

##############################################################################
# We get one SRS file per day. To read this file, we pass the filename into
# the SRS reader. So now `srs_table` contains an astropy table.
srs_table = srs.read_srs(srs_downloaded_files[0])
print(srs_table)

##############################################################################
# We only need the rows which have 'ID' = 'I' or 'IA'.

if 'I' in srs_table['ID'] or 'IA' in srs_table['ID']:
    srs_table = srs_table[np.logical_or(srs_table['ID'] == 'I',
                                        srs_table['ID'] == 'IA')]
else:
    print("Warning : No I or IA entries for this date.")
    srs_table = None

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
                    xycoords=ax.get_transform('heliographic_stonyhurst'))
plt.show()
