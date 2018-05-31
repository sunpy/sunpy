"""
=========================================================
Overplotting Active Region locations on magnetogram Plots
=========================================================

This example shows how to overplot Active Region location on magnetogram plots.
"""

##############################################################################
# Start by importing the necessary modules.
from __future__ import print_function, division

import datetime

import numpy as np
import matplotlib.pyplot as plt

from astropy import units as u
from astropy.coordinates import SkyCoord

import sunpy.map
import sunpy.coordinates
from sunpy.io.special import srs
from sunpy.time import parse_time
from sunpy.net import Fido, attrs as a

##############################################################################
# Let's select a date (yyyy-mm-dd) for which we will be downloading files.

day = parse_time("2017-01-25")

##############################################################################
# We will select the entire day as our timerange.

start_time = day
end_time = day + datetime.timedelta(hours=23, minutes=59, seconds=59)

##############################################################################
# Send the search query.

results = Fido.search(a.Time(start_time, end_time),
                      a.Instrument('HMI') & a.vso.Physobs("LOS_magnetic_field"),
                      a.vso.Sample(60 * u.second))

##############################################################################
# We will only download the first file for the day. For that we use fido
# indexing on the search results which will return the first file for the day.

result = results[0, 0]

##############################################################################
# Download the file. The `fetch` method returns a list of filenames. As we
# used indexing to get the first file of the day, the list contains one
# filename.

file_name = Fido.fetch(result)

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
# Now to make the plot.
# Create the magnetogram plot using the FITS file.

smap = sunpy.map.Map(file_name)

ax = plt.subplot(projection=smap)

smap.plot(vmin=-120, vmax=120)

smap.draw_limb()

ax.set_autoscale_on(False)

##############################################################################
# We make a SkyCoord object and plot the active points on the map.
# Add the numbers as labels for each point.

if len(lats) > 0:
    c = SkyCoord(lngs, lats, frame="heliographic_stonyhurst")
    ax.plot_coord(c, 'o')

    for i, num in enumerate(numbers):
        ax.annotate(num, (lngs[i].value, lats[i].value),
                    xycoords=ax.get_transform('heliographic_stonyhurst'))

##############################################################################
# Now we display the combined plot.

plt.show()
