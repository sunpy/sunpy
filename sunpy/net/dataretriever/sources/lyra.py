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
    ...                       a.Instrument.lyra, a.Level(2))  #doctest: +REMOTE_DATA
    >>> results  #doctest: +REMOTE_DATA
    <sunpy.net.fido_factory.UnifiedResponse object at ...>
    Results from 1 Provider:
    <BLANKLINE>
    2 Results from the LYRAClient:
         Start Time     Source Provider  Physobs   Instrument Level
    ------------------- ------ -------- ---------- ---------- -----
    2016-01-01 00:00:00 Proba2      esa irradiance       lyra     2
    2016-01-02 00:00:00 Proba2      esa irradiance       lyra     2
    <BLANKLINE>
    <BLANKLINE>

    """
    baseurl = (r'http://proba2.oma.be/lyra/data/bsd/%Y/%m/%d/'
               r'lyra_%Y%m%d-000000_lev(\w){1}_std.fits')
    extractor = ('http://proba2.oma.be/lyra/data/bsd/{4d}/{2d}/{2d}'
                 '/lyra_{}-000000_lev{Level:d}_std.fits')

    def _makeimap(self):
        """
        Helper Function:used to hold information about source.
        """
        self.map_['source'] = 'Proba2'
        self.map_['instrument'] = 'lyra'
        self.map_['physobs'] = 'irradiance'
        self.map_['provider'] = 'esa'

    @classmethod
    def _can_handle_query(cls, *query):
        """
        Answers whether client can service the query.

        Parameters
        ----------
        query : list of query objects

        Returns
        -------
        boolean
            answer as to whether client can service the query
        """
        chkattr = ['Time', 'Instrument', 'Level']
        chklist = [x.__class__.__name__ in chkattr for x in query]
        for x in query:
            if x.__class__.__name__ == 'Instrument' and x.value.lower() == 'lyra':
                return all(chklist)
        return False

    @classmethod
    def register_values(cls):
        from sunpy.net import attrs
        adict = {attrs.Instrument: [('LYRA',
                                     'Lyman Alpha Radiometer is the solar UV radiometer on board Proba-2.')],
                 attrs.Level: [('1', 'LYRA: Metadata and uncalibrated data daily fits.'),
                               ('2', 'LYRA: Calibrated data, provided as daily fits.'),
                               ('3', 'LYRA: Same as level 2 but the calibrated data is averaged over 1 min.')]}
        return adict
