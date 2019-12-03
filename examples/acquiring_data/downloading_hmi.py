"""
===========================================
Downloading and plotting an HMI magnetogram
===========================================

This example shows how to download a HMI magnetogram data with Fido and make a plot.
"""
import astropy.units as u
import matplotlib.pyplot as plt

import sunpy.map
from sunpy.net import Fido, attrs as a

###############################################################################
# To download the required data, we use
# `Fido <sunpy.net.fido_factory.UnifiedDownloaderFactory>`, a downloader client,
# to query the Joint Science Operations Center, or JSOC, where HMI data are stored.
# First define the search variables, a timerange, 
# a [data series](http://jsoc.stanford.edu/JsocSeries_DataProducts_map.html),
# keywords, and your e-mail address (to notify you when the download is complete).
# See the JSOC e-mail address registration page
# [here](http://jsoc.stanford.edu/ajax/register_email.html).

result = Fido.search(a.Time('2014/11/20 00:00:00', '2014/11/20 00:04:00'),
                     a.jsoc.Series("hmi.M_720s"),
                     a.jsoc.Keys(["T_REC, CROTA2"]),
                     a.jsoc.Notify("sunpy@sunpy.org"))

###############################################################################
# Now we can see what results we obtained from our search.
# Notice we have two files.
print(result)

###############################################################################
# The following shows how to download the results. If we
# don't provide a path it will download the file into the sunpy data directory.
# The output provides the path of the downloaded files.
# The result can be from several data clients, so we have to index the client first
# and then index the file.

# Slice the first record returned by the first client.
downloaded_file = Fido.fetch(result[0])
print(downloaded_file)

###############################################################################
# Now load it into a map and plot it.
# We see that solar North is pointed down instead of up in this image, which is
# indicated by the coordinates (that range from positive to negative, rather
# than negative to positive).
hmi_map = sunpy.map.Map(downloaded_file[0])
fig = plt.figure()
hmi_map.plot()
plt.show()

###############################################################################
# Rotate the image such that solar North is pointed up.
# The order specifies the type of interpolation; in this case, 3 refers to 
# bi-cubic.
hmi_rotated = hmi_map.rotate(order=3)
hmi_rotated.plot()
plt.show()
