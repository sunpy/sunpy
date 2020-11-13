"""
==========================================================
Overplotting SRS active region locations on a magnetograms
==========================================================

How to find and plot the location of an active region on an HMI magnetogram.
"""
import matplotlib.pyplot as plt
import numpy as np

import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.time import TimeDelta

import sunpy.coordinates
import sunpy.data.sample
import sunpy.map
from sunpy.io.special import srs
from sunpy.net import Fido
from sunpy.net import attrs as a
from sunpy.time import parse_time

##############################################################################
# For this example, we will start with the sample data. We need an HMI
# file and use it to create a map.

smap = sunpy.map.Map(sunpy.data.sample.HMI_LOS_IMAGE)

##############################################################################
# Then using the observation time of the sample image we will download the corresponding
# NOAA SWPC solar region summary. These are published once a day, so we will define an
# ``end_time`` that ends just before the next day, otherwise you will download two files.

start_time = parse_time(smap.date)
end_time = start_time + TimeDelta(23*u.hour + 59*u.minute + 59*u.second)

##############################################################################
#  Here we use `sunpy.net.Fido` to acquire the data.

srs_results = Fido.search(a.Time(start_time, end_time), a.Instrument.srs_table)
srs_downloaded_files = Fido.fetch(srs_results)

##############################################################################
# To read this file, we pass the filename into the SRS reader. So now
# `srs_table` contains an astropy table.
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
