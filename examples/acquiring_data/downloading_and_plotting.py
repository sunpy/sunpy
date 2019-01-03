"""
=============================
Downloading and plotting data
=============================
"""
###############################################################################
# Import the required modules.
from sunpy.net import Fido, attrs
import sunpy.map
import matplotlib.pyplot as plt

###############################################################################
# Define the timerange and instrument.
timerange = attrs.Time('1996/05/24 11:00', '1996/05/24 12:00')
instrument = attrs.Instrument('C3')

###############################################################################
# Search for data that is from the requested instrument and in the requested
# timerange.
result = Fido.search(timerange, instrument)
print(result)

###############################################################################
# Download the files. "downloaded_files" contains a list of the location of
# each of the locally downloaded files.
downloaded_files = Fido.fetch(result)
print(downloaded_files)

###############################################################################
# Create a SunPy Map object from the first downloaded file.
lascomap = sunpy.map.Map(downloaded_files[0])

###############################################################################
# Plot the Map.
fig, ax = plt.subplots()
lascomap.plot(axes=ax)
plt.show()
