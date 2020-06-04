#  Author: Rishabh Sharma <rishabh.sharma.gunner@gmail.com>
#  This Module was developed under funding provided by
#  Google Summer of Code 2014

import astropy.units as u
from astropy.time import TimeDelta

from sunpy.net.dataretriever import GenericClient
from sunpy.time import TimeRange
from sunpy.util.scraper import Scraper

__all__ = ['NoRHClient']

BASEURL = 'ftp://solar-pub.nao.ac.jp/pub/nsro/norh/data/tcx/%Y/%m/{freq}%y%m%d'


class NoRHClient(GenericClient):
    """
    Provides access to the Nobeyama RadioHeliograph (NoRH) averaged correlation
    time series data.

    Uses this `ftp archive <ftp://solar-pub.nao.ac.jp/pub/nsro/norh/data/tcx/>`__
    hosted by the `NoRH Science Center <https://solar.nro.nao.ac.jp/norh/doc/manuale/node1.html>`__.

    Queries to NoRH should specify either 17GHz or 34GHz as a Wavelength.

    Examples
    --------
    >>> import astropy.units as u
    >>> from sunpy.net import Fido, attrs as a
    >>> results = Fido.search(a.Time("2016/1/1", "2016/1/2"),
    ...                       a.Instrument.norh, a.Wavelength(17*u.GHz))  #doctest: +REMOTE_DATA
    >>> results  #doctest: +REMOTE_DATA
    <sunpy.net.fido_factory.UnifiedResponse object at ...>
    Results from 1 Provider:
    <BLANKLINE>
    2 Results from the NoRHClient:
         Start Time           End Time      Source Instrument   Wavelength
    ------------------- ------------------- ------ ---------- --------------
    2016-01-01 00:00:00 2016-01-02 00:00:00   NAOJ       NORH 17000000.0 kHz
    2016-01-02 00:00:00 2016-01-03 00:00:00   NAOJ       NORH 17000000.0 kHz
    <BLANKLINE>
    <BLANKLINE>

    """

    def _get_url_for_timerange(self, timerange, **kwargs):
        """
        Returns list of URLS corresponding to value of input timerange.

        Parameters
        ----------
        timerange: `sunpy.time.TimeRange`
            time range for which data is to be downloaded.

        Returns
        -------
        urls : list
            list of URLs corresponding to the requested time range
        """

        # We allow queries with no Wavelength but error here so that the query
        # does not get passed to VSO and spit out garbage.
        if 'wavelength' not in kwargs.keys() or not kwargs['wavelength']:
            raise ValueError("Queries to NORH should specify either 17GHz or 34GHz as a Wavelength."
                             "see https://solar.nro.nao.ac.jp/norh/doc/manuale/node65.html")
        else:
            wavelength = kwargs['wavelength']

        # If wavelength is a single value GenericClient will have made it a
        # Quantity in the kwargs.
        if not isinstance(wavelength, u.Quantity):
            raise ValueError(f"Wavelength to NORH must be one value not {wavelength}.")

        wavelength = wavelength.to(u.GHz, equivalencies=u.spectral())
        if wavelength == 34 * u.GHz:
            freq = 'tcz'
        elif wavelength == 17 * u.GHz:
            freq = 'tca'
        else:
            raise ValueError("NORH Data can be downloaded for 17GHz or 34GHz,"
                             " see https://solar.nro.nao.ac.jp/norh/doc/manuale/node65.html")

        # If start of time range is before 00:00, converted to such, so
        # files of the requested time ranger are included.
        # This is done because the archive contains daily files.
        if timerange.start.strftime('%M-%S') != '00-00':
            timerange = TimeRange(timerange.start.strftime('%Y-%m-%d'),
                                  timerange.end)

        norh = Scraper(BASEURL, freq=freq)
        # TODO: warn user that some files may have not been listed, like for example:
        #       tca160504_224657 on ftp://solar-pub.nao.ac.jp/pub/nsro/norh/data/tcx/2016/05/
        #       as it doesn't follow pattern.

        return norh.filelist(timerange)

    def _get_time_for_url(self, urls):
        freq = urls[0].split('/')[-1][0:3]  # extract the frequency label
        crawler = Scraper(BASEURL, freq=freq)
        times = list()
        for url in urls:
            t0 = crawler._extractDateURL(url)
            # hard coded full day as that's the normal.
            times.append(TimeRange(t0, t0 + TimeDelta(1*u.day)))
        return times

    def _makeimap(self):
        """
        Helper Function used to hold information about source.
        """
        self.map_['source'] = 'NAOJ'
        self.map_['provider'] = 'NRO'
        self.map_['instrument'] = 'NORH'
        self.map_['physobs'] = ''

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
        # Import here to prevent circular imports
        from sunpy.net import attrs as a

        required = {a.Time, a.Instrument}
        optional = {a.Wavelength}
        if not cls.check_attr_types_in_query(query, required, optional):
            return False

        # if we get this far we have either Instrument and Time
        # or Instrument, Time and Wavelength
        for x in query:
            if isinstance(x, a.Instrument) and x.value.lower() == 'norh':
                return True

        return False

    @classmethod
    def register_values(cls):
        from sunpy.net import attrs
        adict = {attrs.Instrument: [('NORH',
                                     ('Nobeyama Radio Heliograph is an imaging radio telescope at 17 '
                                      'or 34GHz located at the Nobeyama Solar Radio Observatory.'))]}
        return adict
