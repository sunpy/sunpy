import datetime
import numpy as np

import astropy.units as u

from sunpy.net.dataretriever.client import GenericClient
from sunpy.util.scraper import Scraper
from sunpy.net import attrs as a

__all__ = ['KanzelhoheClient']


class KanzelhoheClient(GenericClient):
    """
    Returns a list of URLS to Kanzelhohe H-alpha files corresponding to value of input timerange.
    URL source: `http://cesar.kso.ac.at/`.
    The earliest data for H-alpha 2k - 20-Jul-2000
                          Ca-II k - 31-Jul-2010
                          Continuum - 7-Jan-2011

    Parameters
    ----------
    timerange: `~sunpy.time.TimeRange`
        time range for which data is to be downloaded.
        Example value - TimeRange('2015-12-30 00:00:00','2015-12-31 00:01:00')
    Instrument: Fixed argument = 'kanzelhohe'
    Wavelength: Fixed argument = astropy.units.quantity.Quantity
                The physical value of wavelength will belong to any of [5460, 6563, 32768]
                and units will be Angstroms.

    Returns
    -------
    urls: `list`
    list of urls corresponding to requested time range.
    Examples
    --------
    >>> from sunpy.net import Fido, attrs as a
    >>> timerange = a.Time('2015/12/28 00:00:00','2015/12/28 09:03:00')
    >>> results = Fido.search(timerange, a.Instrument('kanzelhohe'),
    ...                       a.Wavelength(6563*u.AA))   #doctest: +REMOTE_DATA
    >>> print(results)  #doctest: +REMOTE_DATA
    Results from 1 Provider:
    <BLANKLINE>
    1 Results from the KanzelhoheClient:
         Start Time           End Time      ...   Instrument      Wavelength
           str19               str19        ...     str14           str15
    ------------------- ------------------- ... -------------- ---------------
    2015-12-28 00:00:00 2015-12-28 09:03:00 ... Kanzelhohe HA2 6563.0 Angstrom
    <BLANKLINE>
    <BLANKLINE>
    >>> response = Fido.fetch(results)  #doctest: +SKIP
    """

    def _get_url_for_timerange(self, timerange, **kwargs):
        """
        Return list of URLS corresponding to value of input timerange.

        Parameters
        ----------
        timerange: `~sunpy.time.TimeRange`
            time range for which data is to be downloaded.

        Returns
        -------
        urls : `list`
            list of URLs corresponding to the requested time range
        """
        wave = kwargs['wavelength']
        wave_float = (wave.min if isinstance(wave, a.Wavelength) else wave.min()).value
        table = {6563: {'datatype': ['halpha2k/recent', 'halph_fr'],
                        'start_time': datetime.datetime(2000, 7, 20, 7, 45, 46)},
                 32768: {'datatype': ['caiia', 'caiik_fi'],
                         'start_time': datetime.datetime(2000, 7, 20, 7, 45, 46)},
                 5460: {'datatype': ['phokada', 'bband_fi'],
                        'start_time': datetime.datetime(2000, 7, 20, 7, 45, 46)}}
        # Checking if value is close enough to a wavelength value.
        # Converting from one unit to other introduces precision errors.
        wavelengths = table.keys()
        wave_list = list(filter(lambda x: np.isclose(x, wave_float), wavelengths))
        if(len(wave_list) == 0):
            raise ValueError("no value is enough closer to allowed wavelengths")
        wave = wave_list[0]
        start_date = table[wave]['start_time']
        if timerange.start < start_date:
            msg = 'Earliest date for which Kanzelhohe data is available is {:%Y-%m-%d}'
            raise ValueError(msg.format(start_date))
        prefix = "http://cesar.kso.ac.at/{datatype}/%Y/"
        suffix = ""
        if wave != 6563:
            suffix = "%Y%m%d/processed/"
        url_pattern = prefix + suffix + "kanz_{datatype1}_%Y%m%d_%H%M%S.fts.gz"
        dtypes = table[wave]['datatype']
        scraper = Scraper(url_pattern, datatype=dtypes[0], datatype1=dtypes[1])
        if not timerange:
            return []
        result = scraper.filelist(timerange)
        return result

    def _makeimap(self):
        """
        Helper Function:used to hold information about source.
        """
        self.map_['source'] = 'Global Halpha Network'  # TODO: check with kanzelhohe
        self.map_['instrument'] = 'Kanzelhohe HA2'
        self.map_['phyobs'] = 'irradiance'
        self.map_['provider'] = 'Kanzelhohe'

    @classmethod
    def _can_handle_query(cls, *query):
        """
        Answers whether client can service the query.

        Parameters
        ----------
        query : list of query objects

        Returns
        -------
        boolean: answer as to whether client can service the query
        """
        chk_var = 0
        for x in query:
            if x.__class__.__name__ == 'Instrument' and isinstance(
                    x.value, str) and x.value.lower() == 'kanzelhohe':
                chk_var += 1
            if x.__class__.__name__ == 'Wavelength':
                chk_var += 1
        return chk_var == 2
