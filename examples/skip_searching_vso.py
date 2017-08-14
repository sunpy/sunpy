"""
=================
Searching the VSO
=================

A simple example showing how to download data from the VSO with Fido.
"""

###############################################################################
# Fido is the primary interface to search for and download data and
# will search the VSO if appropriate. First import it and the search
# attributes.
from __future__ import print_function, division

import astropy.units as u

from sunpy.net import Fido, attrs as a

###############################################################################
# We could ask for all SOHO/EIT data between January 1st and 2nd, 2001.

attrs_time = a.Time('2005/01/01 00:10', '2005/01/01 00:15')
result = Fido.search(attrs_time, a.Instrument('eit'))

###############################################################################
# Let's inspect the result

print(result)

###############################################################################
# Now lets download this query. If we don't provide a path it will download the
# file into the sunpy data directory.

downloaded_files = Fido.fetch(result)

###############################################################################
# You can check where the file was downloaded to.

print(downloaded_files)

###############################################################################
# More complicated queries can be constructed by using relational operators.
# For example, say we are interested in both eit and mdi data.

result = Fido.search(a.Time('2012/3/4', '2012/3/6'),
                     a.Instrument('aia'),
                     a.Wavelength(171*u.angstrom) | a.Wavelength(94*u.angstrom))
print(result)
