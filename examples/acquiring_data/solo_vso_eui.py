"""
======================================
Downloading EUI Data with VSO and SOAR
======================================

This example demonstrates querying for Solar Orbiter EUI data with
both the SOAR client and VSO client.
"""

# sphinx_gallery_tags = ["Acquiring Data", "Solar Orbiter", "VSO", "SOAR"]

import matplotlib.pyplot as plt

import astropy.units as u

import sunpy.map
from sunpy.net import Fido
from sunpy.net import attrs as a

###############################################################################
# SunPy's Fido supports querying many different sources of data simultaneously.
# This query will search both the Solar Orbiter Archive (SOAR) and the VSO for Solar Orbiter data products.
#
# Some Solar Orbiter products are available through both services, so this
# query may return overlapping results from the two clients.
# This example demonstrates a query which returns (the same) data from both sources and how to select between them.
# Using `a.Provider <sunpy.net.attrs.Provider>` it is possible to restrict a Fido search to only one provider.
combined_results = Fido.search(a.Time("2022-03-01", "2022-03-01 00:10"),
                               a.Instrument("EUI"),
                               a.Wavelength(17.4*u.nm),
                               a.Level(2))

print(combined_results)

################################################################################
# We can select one of the two clients like this:
# VSO:
print(combined_results["vso"])

################################################################################
# SOAR:
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
