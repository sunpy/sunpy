
"""
'attrs' are parameters which can be composed together to specify searches to
`sunpy.net.Fido`. They can be combined
by using logical and (``&``) and logical or (``|``) operations to construct
very complex queries.

For example you could combine two instruments using or (``|``) with a time
specification and a sample cadence using:

.. code-block:: python

    >>> import astropy.units as u
    >>> from sunpy.net import Fido, attrs as a
    >>> Fido.search(a.Time("2011/01/01", "2011/01/02") & (a.Instrument.aia | a.Instrument.hmi) & a.Sample(1*u.day))  # doctest: +REMOTE_DATA

In addition to the core attrs defined here, other sunpy clients also provide
attrs specific to them, under:

* `a.vso <sunpy.net.vso.attrs>`
* `a.jsoc <sunpy.net.jsoc.attrs>`
* `a.goes <sunpy.net.dataretriever.attrs.goes>`

"""
from ._attrs import Detector, Instrument, Level, Physobs, Resolution, Sample, Time, Wavelength

# Trick the docs into thinking these attrs are defined in here.
for _a in (Time, Instrument, Wavelength, Level, Sample, Detector, Resolution, Physobs):
    _a.__module__ = __name__

__all__ = ['Time', 'Instrument', 'Wavelength', 'Level',
           'Sample', 'Detector', 'Resolution', 'Physobs']
