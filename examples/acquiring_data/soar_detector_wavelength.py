"""
=========================================================================
Searching for Solar Orbiter data using Wavelength and Detector attributes
=========================================================================

This example demonstrates how to search and download Solar Orbiter data using ``sunpy.net.Fido``.
To do this, we can build a query based on several attributes.
Here, we will build a search for METIS data from the UVD (Ultra Violet Detector) detector for a specific wavelength.
"""

import astropy.units as u

import sunpy.net.attrs as a
from sunpy.net import Fido

###############################################################################
# We shall start with constructing a search query with wavelength and detector.

instrument = a.Instrument("METIS")
time = a.Time("2023-02-01 01:00", "2023-02-01 05:00")
level = a.Level(2)
detector = a.Detector("UVD")
wavelength = a.Wavelength(121.6 * u.AA)

###############################################################################
# Now do the search.

result = Fido.search(instrument & time & level & detector & wavelength)
result

###############################################################################
# To then download the data, you would then use Fido.fetch(result), which will download the data locally.
