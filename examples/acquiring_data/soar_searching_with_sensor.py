"""
===============================================================
Searching for EPD/EPT data with the ``a.soar.Sensor`` attribute
===============================================================

This example demonstrates how to search and download Solar Orbiter data for a specific instrument sensor.
Here, we will build a query for EPD data, specifically from the EPT sensor using `a.soar.Sensor <sunpy.net.soar.attrs.Sensor>`.
"""

import sunpy.net.attrs as a
from sunpy.net import Fido

###############################################################################
# Importing sunpy.net.soar registers the client with sunpy Fido

import sunpy.net.soar  # NOQA: F401 isort:skip

###############################################################################
# We shall construct a search query with instrument, time, level, and sensor attributes.

result = Fido.search(
    a.Instrument("EPD"),
    a.Time("2023-02-01", "2023-02-03"),
    a.Level(2),
    a.soar.Sensor("EPT"),
)

###############################################################################
# And we can see the result as follows:

result

###############################################################################
# To then download the data, you would then use ``Fido.fetch(result)``, which will download the data locally.
