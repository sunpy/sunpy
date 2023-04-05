# Author: Rishabh Sharma <rishabh.sharma.gunner@gmail.com>
# This module was developed under funding provided by
# Google Summer of Code 2014

from sunpy.net.dataretriever.client import GenericClient

__all__ = ['LYRAClient']


class LYRAClient(GenericClient):
    """
    Provides access to the LYRA/Proba2 data archive.

    Hosted by the `PROBA2 Science Center <http://proba2.oma.be/lyra/data/bsd/>`__.

    Examples
    --------
    >>> from sunpy.net import Fido, attrs as a
    >>> results = Fido.search(a.Time("2016/1/1", "2016/1/2"),
    ...                       a.Instrument.lyra)  #doctest: +REMOTE_DATA
    >>> results  #doctest: +REMOTE_DATA
    <sunpy.net.fido_factory.UnifiedResponse object at ...>
    Results from 1 Provider:
    <BLANKLINE>
    4 Results from the LYRAClient:
    Source: http://proba2.oma.be/lyra/data/bsd
    <BLANKLINE>
           Start Time               End Time        Source ... Extent Type   Size
                                                           ...              Mibyte
    ----------------------- ----------------------- ------ ... ----------- --------
    2016-01-01 09:41:00.000 2016-01-01 10:40:00.000 PROBA2 ...         N/A  2328.75
    2016-01-01 09:41:00.000 2016-01-01 10:40:00.000 PROBA2 ...         N/A 419.0625
    2016-01-01 09:41:00.000 2016-01-01 10:40:00.000 PROBA2 ...         N/A  30.9375
    <BLANKLINE>
    """
    baseurl = (r'http://proba2.oma.be/lyra/data/bsd/%Y/%m/%d/'
               r'lyra_(\d){8}-000000_lev(\d){1}_std.fits')
    pattern = '{}/bsd/{year:4d}/{month:2d}/{day:2d}/{}_lev{Level:1d}_std.fits'

    @property
    def info_url(self):
        return 'http://proba2.oma.be/lyra/data/bsd'

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
