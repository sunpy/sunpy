# -*- coding: utf-8 -*-
"""
======================================
Downloading and plotting LASCO C3 data
======================================

How to download SOHO/LASCO C3 data with Fido and make a plot.
"""
import matplotlib.pyplot as plt

from sunpy.net import Fido, attrs
import sunpy.map
from sunpy.io.file_tools import read_file

###############################################################################
# In order to download the required data, we use
# `Fido <sunpy.net.fido_factory.UnifiedDownloaderFactory>`, a downloader client.
# We define two search variables:
# a timerange and the instrument.
timerange = attrs.Time('1998/05/24 11:00', '1998/05/24 11:20')
instrument = attrs.Instrument('LASCO')
detector = attrs.Detector('C3')
result = Fido.search(timerange, instrument)

###############################################################################
# Let's inspect the result
print(result)

###############################################################################
# The following shows how to download the results. If we
# don't provide a path it will download the file into the sunpy data directory.
# The output provides the path of the downloaded files.
downloaded_files = Fido.fetch(result[0])
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
lascomap = sunpy.map.Map(data, header)
fig = plt.figure()
lascomap.plot()
plt.show()
