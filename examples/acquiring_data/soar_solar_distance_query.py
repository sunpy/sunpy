"""
===============================================================
Searching for Solar Orbiter data using Solar Distance attribute
===============================================================

This example demonstrates how to search and download Solar Orbiter data using ``sunpy.net.Fido``.
To do this, we can build a query based on several attributes.

The ``Distance`` attribute allows us to specify a range of distances from the Sun in astronomical units (AU).

"""

# sphinx_gallery_tags = ["Acquiring Data", "SOAR", "Solar Orbiter"]

import astropy.units as u

from sunpy.net import Fido
from sunpy.net import attrs as a

###############################################################################
# We shall start with constructing a search query with instrument, level, detector, and distance.

instrument = a.Instrument("EUI")
time = a.Time("2022-10-29 05:00:00", "2022-10-29 06:00:00")
level = a.Level(2)
detector = a.Detector("HRI_EUV")
distance = a.soar.Distance(0.45 * u.AU, 0.46 * u.AU)

###############################################################################
# Now do the search without time attribute.

result = Fido.search(instrument & level & detector & distance)
result

###############################################################################
# Now do the search with time attribute.
result = Fido.search(instrument & level & detector & distance & time)
result

###############################################################################
# To then download the data, you would use Fido.fetch(result), which will download the data locally.
