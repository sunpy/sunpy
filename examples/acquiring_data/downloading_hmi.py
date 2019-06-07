"""
===========================================
Downloading and plotting an HMI magnetogram
===========================================

How to download an HMI magnetogram data with Fido and make a plot.
"""
import astropy.units as u
import matplotlib.pyplot as plt

import sunpy.map
from sunpy.net import Fido, attrs as a

###############################################################################
# To download the required data, we use
# `Fido <sunpy.net.fido_factory.UnifiedDownloaderFactory>`, a downloader client.
# First define the search variables, a timerange, the instrument,
# the observation type,
# and a cadence of images spaced every 720 seconds.
result = Fido.search(a.Time('2011/11/09 17:40:00', '2011/11/09 17:55:00'),
                     a.Instrument('hmi'),
                     a.Sample(720*u.s),
                     a.vso.Physobs('LOS_magnetic_field'))

###############################################################################
# Now we can see what results we obtained from our search.
# Notice we have two files.
print(result)

###############################################################################
# The following shows how to download the results. If we
# don't provide a path it will download the file into the sunpy data directory.
# The output provides the path of the downloaded files.
# The result can be from several data clients, so we have to index the client first
# client and then index the file.

# Slice the first record returned by the first client.
downloaded_file = Fido.fetch(result[0, 0])
print(downloaded_file)

###############################################################################
# Now load it into a map and plot it
hmi_map = sunpy.map.Map(downloaded_file[0])
fig = plt.figure()
hmi_map.plot()
plt.show()
