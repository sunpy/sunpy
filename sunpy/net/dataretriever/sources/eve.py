# Author: Rishabh Sharma <rishabh.sharma.gunner@gmail.com>
# This module was developed under funding by
# Google Summer of Code 2014

from sunpy.net.dataretriever import GenericClient

__all__ = ['EVEClient']


class EVEClient(GenericClient):
    """
    Provides access to Level 0C Extreme ultraviolet Variability Experiment (EVE) data.

    To use this client you must request Level 0 data.
    It is hosted by `LASP <http://lasp.colorado.edu/home/eve/data/data-access/>`__.

    Examples
    --------

    >>> from sunpy.net import Fido, attrs as a
    >>> results = Fido.search(a.Time("2016/1/1", "2016/1/2"),
    ...                       a.Instrument.eve, a.Level.zero)  #doctest: +REMOTE_DATA
    >>> results  #doctest: +REMOTE_DATA
    <sunpy.net.fido_factory.UnifiedResponse object at ...>
    Results from 1 Provider:
    <BLANKLINE>
    2 Results from the EVEClient:
           Start Time               End Time        Instrument ... Provider Level
    ----------------------- ----------------------- ---------- ... -------- -----
    2016-01-01 00:00:00.000 2016-01-01 23:59:59.999        EVE ...     LASP     0
    2016-01-02 00:00:00.000 2016-01-02 23:59:59.999        EVE ...     LASP     0
    <BLANKLINE>
    <BLANKLINE>

    """
    baseurl = (r'http://lasp.colorado.edu/eve/data_access/evewebdata/quicklook/'
               r'L0CS/SpWx/%Y/%Y%m%d_EVE_L0CS_DIODES_1m.txt')
    pattern = '{}/SpWx/{:4d}/{year:4d}{month:2d}{day:2d}_EVE_L{Level:1d}{}'

    @classmethod
    def register_values(cls):
        from sunpy.net import attrs
        adict = {attrs.Instrument: [('EVE', 'Extreme ultraviolet Variability Experiment, which is part of the NASA Solar Dynamics Observatory mission.')],
                 attrs.Physobs: [('irradiance', 'the flux of radiant energy per unit area.')],
                 attrs.Source: [('SDO', 'The Solar Dynamics Observatory.')],
                 attrs.Provider: [('LASP', 'The Laboratory for Atmospheric and Space Physics.')],
                 attrs.Level: [('0', 'EVE: The specific EVE client can only return Level 0C data. Any other number will use the VSO Client.')]}
        return adict
