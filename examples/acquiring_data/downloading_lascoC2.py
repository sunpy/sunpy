"""
======================================
Downloading and plotting LASCO C2 data
======================================

How to download SOHO/LASCO C2 data with Fido and make a plot.
"""
import matplotlib.pyplot as plt

import sunpy.map
from sunpy.net import Fido
from sunpy.net import attrs as a

###############################################################################
# In order to download the required data, we use `sunpy.net.Fido`, a downloader client.
# We define two search variables: a timerange and the instrument.

timerange = a.Time('2002/05/24 11:06', '2002/05/24 11:07')
instrument = a.Instrument.lasco
detector = a.Detector.c2
result = Fido.search(timerange, instrument, detector)

###############################################################################
# Let's inspect the result.
print(result)

###############################################################################
# The following shows how to download the results. If we
# don't provide a path it will download the file into the sunpy data directory.
# The output provides the path of the downloaded files.

downloaded_files = Fido.fetch(result)
print(downloaded_files)

###############################################################################
# We can then load a downloaded file into a sunpy map and plot it.

lascomap = sunpy.map.Map(downloaded_files[0])
fig = plt.figure()
lascomap.plot()

plt.show()
