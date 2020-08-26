"""
======================================
Searching and downloading from the VSO
======================================

How to download data from the VSO with Fido.
"""
import astropy.units as u

from sunpy.net import Fido
from sunpy.net import attrs as a

###############################################################################
# `sunpy.net.Fido` is the primary
# interface to search for and download data and
# will search the VSO when appropriate. The following example searches for all
# SOHO/EIT images between the times defined below by defining
# a timerange (`~sunpy.net.attrs.Time`) and
# the instrument (`~sunpy.net.attrs.Instrument`).
attrs_time = a.Time('2005/01/01 00:10', '2005/01/01 00:15')
result = Fido.search(attrs_time, a.Instrument.eit)

###############################################################################
# Let's inspect the results.
print(result)

###############################################################################
# The following shows how to download the results. If we
# don't provide a path it will download the file into the sunpy data directory.
# The output provides the path of the downloaded files.
downloaded_files = Fido.fetch(result)
print(downloaded_files)

###############################################################################
# More complicated queries can be constructed by using relational operators.
# For example, it is possible to query two wavelengths at the same time with
# the OR operator (|).
result = Fido.search(a.Time('2020/03/04 00:00', '2020/03/04 00:02'),
                     a.Instrument.aia,
                     a.Wavelength(171*u.angstrom) | a.Wavelength(94*u.angstrom))
print(result)
