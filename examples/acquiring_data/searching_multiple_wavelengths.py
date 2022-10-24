"""
==============================================
Searching for multiple wavelengths with Fido
==============================================
This example shows how you can search for several wavelengths
of AIA data with Fido.
"""
from astropy import units as u

from sunpy.net import Fido, attr
from sunpy.net import attrs as a

###############################################################################
# Here we are demonstrating how you can search for specific wavelengths of data
# from AIA using `Fido` and the `sunpy.net.attr.or_` function.
# For example, you may only want a single wavelength, say 171 Angstrom.:

aia_search = Fido.search(a.Time("2022-02-20 00:00", "2022-02-20 00:02"),
                         a.Instrument("AIA"),
                         a.Wavelength(171*u.angstrom))

aia_search


###############################################################################
# But say you actually want to search for several wavelengths, rather than just one.
# You could use the "|" operator, or instead you can use the `sunpy.net.attr.or_`
# function.

wavelengths = [94, 131, 171, 193, 211]*u.angstrom
aia_search = Fido.search(a.Time("2022-02-20 00:00", "2022-02-20 00:02"),
                         a.Instrument("AIA"),
                         attr.or_(*[a.Wavelength(wav) for wav in wavelengths]))

aia_search

# This returns several searches for each of the wavelengths, which can be indexed.

aia_search[0]
