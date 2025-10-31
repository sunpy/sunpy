# Author: Rishabh Sharma <rishabh.sharma.gunner@gmail.com>
# This module was developed under funding provided by
# Google Summer of Code 2014

from sunpy.net.dataretriever.client import GenericClient

__all__ = ['LYRAClient']


class LYRAClient(GenericClient):
    """
    Provides access to the LYRA/Proba2 data archive.

    Hosted by the `PROBA2 Science Center <https://proba2.sidc.be/lyra>`__.

    Examples
    --------
    >>> from sunpy.net import Fido, attrs as a
    >>> results = Fido.search(a.Time("2025/1/1", "2025/1/2"),
    ...                       a.Instrument.lyra)  # doctest: +REMOTE_DATA
    >>> results  # doctest: +REMOTE_DATA
    <sunpy.net.fido_factory.UnifiedResponse object at ...>
    Results from 1 Provider:
    <BLANKLINE>
    4 Results from the LYRAClient:
    Source: https://proba2.sidc.be/lyra/data/bsd
    <BLANKLINE>
           Start Time               End Time        Instrument  Physobs   Source Provider Level
    ----------------------- ----------------------- ---------- ---------- ------ -------- -----
    2025-01-01 00:00:00.000 2025-01-01 23:59:59.999       LYRA irradiance PROBA2      ESA     2
    2025-01-01 00:00:00.000 2025-01-01 23:59:59.999       LYRA irradiance PROBA2      ESA     3
    2025-01-02 00:00:00.000 2025-01-02 23:59:59.999       LYRA irradiance PROBA2      ESA     2
    2025-01-02 00:00:00.000 2025-01-02 23:59:59.999       LYRA irradiance PROBA2      ESA     3
    <BLANKLINE>
    <BLANKLINE>
    """
    pattern = 'https://proba2.sidc.be/lyra/data/bsd/{{year:4d}}/{{month:2d}}/{{day:2d}}/{{}}_lev{{Level:1d}}_std.fits'

    @property
    def info_url(self):
        return 'https://proba2.sidc.be/lyra/data/bsd'

    @classmethod
    def register_values(cls):
        from sunpy.net import attrs
        adict = {attrs.Instrument: [('LYRA',
                                     'Lyman Alpha Radiometer is the solar UV radiometer on board Proba-2.')],
                 attrs.Physobs: [('irradiance', 'the flux of radiant energy per unit area.')],
                 attrs.Source: [('PROBA2', 'The PROBA-2 Satellite')],
                 attrs.Provider: [('ESA', 'The European Space Agency.')],
                 attrs.Level: [('1', 'LYRA: Metadata and uncalibrated data daily fits.'),
                               ('2', 'LYRA: Calibrated data, provided as daily fits.'),
                               ('3', 'LYRA: Same as level 2 but the calibrated data is averaged over 1 min.')]}
        return adict
