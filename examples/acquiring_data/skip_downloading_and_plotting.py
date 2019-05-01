"""
======================================
Downloading and plotting LASCO C3 data
======================================

In this example, we show how to acquire and plot SOHO/LASCO C3 data using SunPy.
"""
# Import the required modules.
from sunpy.net import Fido, attrs
import sunpy.map
from sunpy.io.file_tools import read_file

###############################################################################
# In order to download the required data, we will use
# `Fido <sunpy.net.fido_factory.UnifiedDownloaderFactory>`, a downloader client.
# Using the `Fido.search` method, we need to define two variables:
# a timerange (`~sunpy.net.attrs.Time`) and
# a instrument (`~sunpy.net.attrs.Instrument`).

timerange = attrs.Time('1996/05/24 11:00', '1996/05/24 12:00')
instrument = attrs.Instrument('C3')

###############################################################################
# This now must be passed in ``Fido.search`` which will query the
# online services:

result = Fido.search(timerange, instrument)

# ``result`` contains the return from the online search.
# In this case, we have found 4 files that correspond to our search parameters.
print(result)

###############################################################################
# The next step is to download the search results and `Fido.fetch` will be used.
# We will pass in the ``result`` and ``downloaded_files`` will contain
# a list of the location of each of the downloaded files.

downloaded_files = Fido.fetch(result)
print(downloaded_files)

###############################################################################
# Finally we can pass in the first file we downloaded into `sunpy.map.Map`
# to create a SunPy Map object which allows us to easily plot the image.
# However, there is a problem here. The downloaded file lacks the correct meta
# information needed to create a `sunpy.map.Map`.

# We want to open the file and access both the data and the header information.
data, header = read_file(downloaded_files[0])[0]

# Add the missing meta information to the header
header['CUNIT1'] = 'arcsec'
header['CUNIT2'] = 'arcsec'

# Now we can pass in the data and header directly into `sunpy.map.Map`
lascomap = sunpy.map.Map(data, header)
lascomap.peek()
