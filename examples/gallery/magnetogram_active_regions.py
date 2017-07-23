"""
=========================================================
Overplotting Active Region locations on magnetogram Plots
=========================================================

This example shows how to overplot Active Region location on magnetogram plots.
"""

##############################################################################
# Start by importing the necessary modules.

import datetime

import matplotlib.pyplot as plt
import numpy as np
from astropy import units as u
from astropy.coordinates import SkyCoord

from sunpy.io.special import srs
import sunpy.coordinates
from sunpy.net import Fido, attrs as a
import sunpy.map
from sunpy.time import parse_time

##############################################################################
# Let's select a date (yyyy-mm-dd) for which we will be downloading files.

day = parse_time("2017-01-01")

##############################################################################
# We will select a small time range to avoid downloading too many files.

start_time = day + datetime.timedelta(minutes=1)
end_time = day + datetime.timedelta(minutes=2)

##############################################################################
# Send the search query.

results = Fido.search(	a.Time(start_time, end_time),
                       a.Instrument('HMI') & a.vso.Physobs(
                           "LOS_magnetic_field"),
                       a.vso.Sample(60 * u.second))

##############################################################################
# Download the files.

downloaded_files = Fido.fetch(results)

##############################################################################
# We will plot only one file in this example.

file_name = downloaded_files[0]

##############################################################################
# Now to download and read the SRS file.
# Download the SRS file.

srs_results = Fido.search(a.Time(start_time, end_time), a.Instrument('SOON'))
srs_downloaded_files = Fido.fetch(srs_results)

##############################################################################
# We get one SRS file per day. So we pass the filename into the SRS reader. So
# now `srs_table` contains an astropy table.

srs_table = srs.read_srs(srs_downloaded_files[0])
print(srs_table)

##############################################################################
# We only need the rows which have 'ID' = 'I' or 'IA'.

if 'I' in srs_table['ID'] or 'IA' in srs_table['ID']:
    srs_table = srs_table[np.logical_or(
        srs_table['ID'] == 'I', srs_table['ID'] == 'IA')]
else:
    print("Warning : No I or IA entries for this date.")
    srs_table = None

##############################################################################
# Now we extract the latitudes, longitudes and the region numbers.

if srs_table is not None:
    lats = srs_table['Latitude']
    lngs = srs_table['Longitude']
    numbers = srs_table['Number']
else:
    lats = lngs = numbers = None

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

c = SkyCoord(lngs, lats, frame="heliographic_stonyhurst")
ax.plot_coord(c, 'o')

##############################################################################
# Add the numbers as labels for each point.

for i, num in enumerate(numbers):
    ax.annotate(num, (lngs[i].value, lats[i].value),
                xycoords=ax.get_transform('heliographic_stonyhurst'))

##############################################################################
# Now we display the combined plot.

plt.show()
