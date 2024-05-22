"""
============================================
Quick overview of using sunpy-soar with Fido
============================================

This example demonstrates how to search and download Solar Orbiter data using ``sunpy.net.Fido``.
"""

import sunpy.net.attrs as a
from sunpy.net import Fido

#####################################################
# Importing sunpy.net.soar registers the client with sunpy
import sunpy.net.soar  # NOQA: F401

#####################################################
# We shall start with constructing a search query.

instrument = a.Instrument("EUI")
time = a.Time("2021-02-01", "2021-02-02")
level = a.Level(1)
product = a.soar.Product("EUI-FSI174-IMAGE")

#####################################################
# Now do the search.

result = Fido.search(instrument & time & level & product)
result

#####################################################
# Finally we can download the data.
#
# For this example, we will comment out the download part
# as we want to avoid downloading data in the documentation build

# files = Fido.fetch(result)
# print(files)
