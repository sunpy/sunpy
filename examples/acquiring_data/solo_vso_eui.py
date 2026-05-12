"""
======================================
Downloading EUI Data with VSO and SOAR
======================================

This example demonstrates querying for Solar Orbiter EUI data with
both the SOAR client and VSO client.
"""

# sphinx_gallery_tags = ["Acquiring Data", "Solar Orbiter", "VSO", "SOAR"]

import matplotlib.pyplot as plt
import sunpy_soar  # noqa

import astropy.units as u

import sunpy.map
from sunpy.net import Fido
from sunpy.net import attrs as a

###############################################################################
# This query (using no specific SOAR attributes) will query both the
# SOAR and the VSO for matching files.
combined_results = Fido.search(a.Time("2022-03-01", "2022-03-01 00:10"),
                                      a.Instrument("EUI"),
                                      a.Wavelength(17.4*u.nm),
                                      a.Level(2))

print(combined_results)

################################################################################
# We can select one of the two clients like this:
print("VSO Only")
print(combined_results["vso"])

print("SOAR Only")
print(combined_results["soar"])


################################################################################
# We shall download one file from the SOAR
filename = Fido.fetch(combined_results["soar"][0])

################################################################################
# Then construct a map and plot it.
euimap = sunpy.map.Map(filename)

euimap.plot()
euimap.draw_grid()
plt.show()
