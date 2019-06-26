# -*- coding: utf-8 -*-
"""
======================================
Downloading and plotting LASCO C3 data
======================================

How to download SOHO/LASCO C3 data with Fido and make a plot.
"""
import matplotlib.pyplot as plt

from sunpy.net import Fido, attrs as a
import sunpy.map
from sunpy.io.file_tools import read_file

###############################################################################
# In order to download the required data, we use
# `Fido <sunpy.net.fido_factory.UnifiedDownloaderFactory>`, a downloader client.
# We define two search variables:
# a timerange and the instrument.
timerange = a.Time('1998/05/24 11:00', '1998/05/24 11:20')
instrument = a.Instrument('LASCO')
detector = a.Detector('C3')
result = Fido.search(timerange, instrument)

###############################################################################
# Let's inspect the result
print(result)

###############################################################################
# The following shows how to download the results. If we
# don't provide a path it will download the file into the sunpy data directory.
# The output provides the path of the downloaded files.
downloaded_files = Fido.fetch(result)
print(downloaded_files)

###############################################################################
# The downloaded file lacks the correct meta. We want to open the file and
# access both the data and the header information.
data, header = read_file(downloaded_files[0])[0]

# Add the missing meta information to the header
header['CUNIT1'] = 'arcsec'
header['CUNIT2'] = 'arcsec'

###############################################################################
# With this fix we can load it into a map and plot the results.
# Please note that there is no plot displayed below as this example is skipped
# due to timeouts that can occur when you try to download LASCO C3 data.
lascomap = sunpy.map.Map(data, header)
fig = plt.figure()
lascomap.plot()
plt.show()
