
"""
'attrs' are parameters which can be composed together to specify searches to
`sunpy.net.Fido`. They can be combined
by using logical and (``&``) and logical or (``|``) operations to construct
very complex queries.

For example you could combine two instruments using or (``|``) with a time
specification and a sample cadence using::

    >>> import astropy.units as u
    >>> from sunpy.net import Fido, attrs as a
    >>> a.Time("2011/01/01", "2011/01/02") & (a.Instrument.aia | a.Instrument.hmi) & a.Sample(1*u.day))  # doctest: +SKIP

In addition to the core attrs defined here, other sunpy clients also provide
attrs specific to them, under:

* `a.adapt <sunpy.net.dataretriever.attrs.adapt>`
* `a.goes <sunpy.net.dataretriever.attrs.goes>`
* `a.hek <sunpy.net.hek.attrs>`
* `a.helio <sunpy.net.helio.attrs>`
* `a.jsoc <sunpy.net.jsoc.attrs>`
* `a.vso <sunpy.net.vso.attrs>`
"""
from ._attrs import (
    Detector,
    ExtentType,
    Instrument,
    Level,
    Physobs,
    Provider,
    Resolution,
    Sample,
    Source,
    Time,
    Wavelength,
)
from .attr import AttrAnd, AttrOr

# Trick the docs into thinking these attrs are defined in here.
for _a in (Time, Instrument, Wavelength, Level, Sample, Detector, Resolution, Physobs, Source, Provider, ExtentType, AttrAnd, AttrOr):
    _a.__module__ = __name__

__all__ = ['Time', 'Instrument', 'Wavelength', 'Level', 'ExtentType',
           'Sample', 'Detector', 'Resolution', 'Physobs', 'Source', 'Provider', 'AttrAnd', 'AttrOr']
