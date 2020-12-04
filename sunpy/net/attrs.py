
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
    <sunpy.net.fido_factory.UnifiedResponse object at ...>
    Results from 2 Providers:
    <BLANKLINE>
    2 Results from the VSOClient:
       Start Time [1]       End Time [1]    Source ...   Type   Wavelength [2]
                                                   ...             Angstrom
    ------------------- ------------------- ------ ... -------- --------------
    2011-01-01 00:00:00 2011-01-01 00:00:01    SDO ... FULLDISK 171.0 .. 171.0
    2011-01-02 00:00:00 2011-01-02 00:00:01    SDO ... FULLDISK 171.0 .. 171.0
    <BLANKLINE>
    11 Results from the VSOClient:
       Start Time [1]       End Time [1]    Source ...   Type    Wavelength [2]
                                                   ...              Angstrom
    ------------------- ------------------- ------ ... -------- ----------------
    2011-01-01 00:00:00                None    SDO ... SYNOPTIC 6173.0 .. 6173.0
    2011-01-01 00:00:00                None    SDO ... SYNOPTIC 6173.0 .. 6173.0
    2011-01-01 00:00:00                None    SDO ... SYNOPTIC 6173.0 .. 6173.0
    2011-01-01 00:00:25 2011-01-01 00:00:26    SDO ... FULLDISK 6173.0 .. 6174.0
    2011-01-01 00:00:25 2011-01-01 00:00:26    SDO ... FULLDISK 6173.0 .. 6174.0
    2011-01-01 00:00:25 2011-01-01 00:00:26    SDO ... FULLDISK 6173.0 .. 6174.0
    2011-01-01 00:10:10 2011-01-01 00:10:11    SDO ... FULLDISK 6173.0 .. 6174.0
    2011-01-01 00:10:10 2011-01-01 00:10:11    SDO ... FULLDISK 6173.0 .. 6174.0
    2011-01-02 00:00:00                None    SDO ... SYNOPTIC 6173.0 .. 6173.0
    2011-01-02 00:00:00                None    SDO ... SYNOPTIC 6173.0 .. 6173.0
    2011-01-02 00:00:00                None    SDO ... SYNOPTIC 6173.0 .. 6173.0
    <BLANKLINE>
    <BLANKLINE>

In addition to the core attrs defined here, other sunpy clients also provide
attrs specific to them, under:

* `a.vso <sunpy.net.vso.attrs>`
* `a.jsoc <sunpy.net.jsoc.attrs>`
* `a.goes <sunpy.net.dataretriever.attrs.goes>`

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

# Trick the docs into thinking these attrs are defined in here.
for _a in (Time, Instrument, Wavelength, Level, Sample, Detector, Resolution, Physobs, Source, Provider, ExtentType):
    _a.__module__ = __name__

__all__ = ['Time', 'Instrument', 'Wavelength', 'Level', 'ExtentType',
           'Sample', 'Detector', 'Resolution', 'Physobs', 'Source', 'Provider']
