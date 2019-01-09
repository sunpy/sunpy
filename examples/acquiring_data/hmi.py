"""
====================================
Downloading and plotting a HMI image
====================================

In this example, we will demonstrate how to download and plot
line of sight magnetic field data from HMI.
"""
# Start by importing the necessary modules.
import astropy.units as u

import sunpy.map
from sunpy.net import Fido, attrs as a

###############################################################################
# Now we will download some data with `sunpy.net.Fido`.
# A `Fido.search` requires us to specify a `~sunpy.net.attr.Time`,
# `~sunpy.net.attr.Sample`, `~sunpy.net.attr.Instrument`
# and the `~sunpy.net.attr.vso.Physobs`.

# We set a time range from ``2015/11/04 12:00:00`` to ``2015/11/04 12:10:00``
# for HMI ``LOS_magnetic_field`` with the images spaced every 720 seconds.
result = Fido.search(a.Time('2015/11/04 12:00:00', '2015/11/04 12:10:00'),
                     a.Instrument('hmi'),
                     a.Sample(720*u.s),
                     a.vso.Physobs('LOS_magnetic_field'))

###############################################################################
# Now we can see what results we obtained from our search.
print(result)

###############################################################################
# Once we are happy with the results obtained from the search.
# We can download the data with `Fido.fetch` .
# In this case we only want one file so we can index the result.

# Notice we have two files. One is the full disk image we plan to display
# and the other is a synoptic version of said image.
downloaded_files = Fido.fetch(result)

###############################################################################
# The ``downloaded_files`` variables returns a list with the path for the file
# that was downloaded.
print(downloaded_files)

###############################################################################
# Now we are going to create a `sunpy.map.Map` object out of the first result.

hmi_map = sunpy.map.Map(downloaded_files[0])

###############################################################################
# Once we have created our Map, we can plot it quite simply doing:

hmi_map.peek()
