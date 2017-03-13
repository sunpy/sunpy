"""
=========================================================
Overplotting Active Region locations on magnetogram Plots
=========================================================

This example shows how to overplot Active Region location on magnetogram plots.
"""

##############################################################################
# Start by importing the necessary modules.

from sunpy.net import Fido, attrs as a
from astropy import units as u
import sunpy.map
import matplotlib.pyplot as plt
import numpy as np

from sunpy.io.special import srs
import sunpy.coordinates
from astropy.coordinates import SkyCoord

##############################################################################
# Let's select a date (yyyy-mm-dd) for which we will be downloading files.

day = "2017-01-01"

##############################################################################
# We will select a small time range to avoid downloading too many files.

start_time = day + " 00:01:00"
end_time = day + " 00:02:00"

##############################################################################
# Send the search query.

results = Fido.search(a.Time(start_time, end_time),
						a.Instrument('HMI') & a.vso.Physobs("LOS_magnetic_field"),
 						a.vso.Sample(60* u.second))

##############################################################################
# Download the files.

downresp = Fido.fetch(results)
#downresp = ['/home/punya/sunpy/data/hmi_m_45s_2017_01_01_00_02_15_tai_magnetogram.fits']

##############################################################################
# We will plot only one file in this example.

file_name = downresp[0]

##############################################################################
# Now to download and read the SRS file.

##############################################################################
# Download the SRS file.

srs_results = Fido.search(a.Time(start_time, end_time), a.Instrument('SOON'))
#srs_downresp = ['/home/punya/sunpy/data/20170101SRS.txt']
srs_downresp = Fido.fetch(srs_results)

##############################################################################
# We get one SRS file per day. So we pass the filename into the SRS reader. So
# now `srs_table` contains an astropy table.

srs_table = srs.read_srs(srs_downresp[0])
print (srs_table)

##############################################################################
# We only need the rows which have 'ID' = 'I' or 'IA'.

srs_table = srs_table[np.logical_or(srs_table['ID'] == 'I', srs_table['ID'] == 'IA')]

##############################################################################
# Now we extract the latitudes, longitudes and the region numbers.
lats = srs_table['Latitude']
lngs = srs_table['Longitude']
numbers = srs_table['Number']

##############################################################################
# Now we make the plot.

##############################################################################
# Create the magnetogram plot using the FITS file.

smap = sunpy.map.Map(file_name)

im = smap.plot()
ax = plt.gca()

##############################################################################
# We make a SkyCoord object and plot the active points on the map.

c = SkyCoord(lngs, lats, frame="heliographic_stonyhurst")
ax.plot_coord(c, 'o')

##############################################################################
# Add the numbers as labels for each point.

for i, num in enumerate(numbers):
	ax.annotate(num, (lngs[i].value, lats[i].value), xycoords=ax.get_transform('heliographic_stonyhurst'))

##############################################################################
# Now we display the combined plot.

smap.plot(vmin=-120, vmax=120)
smap.draw_limb()
plt.show()
