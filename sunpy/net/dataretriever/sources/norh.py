#  Author: Rishabh Sharma <rishabh.sharma.gunner@gmail.com>
#  This Module was developed under funding provided by
#  Google Summer of Code 2014

import datetime
import astropy.units as u

from sunpy.time import TimeRange
from sunpy.util.scraper import Scraper

from sunpy.net import attrs as a
from ..client import GenericClient

__all__ = ['NoRHClient']

BASEURL = 'ftp://solar-pub.nao.ac.jp/pub/nsro/norh/data/tcx/%Y/%m/{freq}%y%m%d'


class NoRHClient(GenericClient):

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
                             "see http://solar.nro.nao.ac.jp/norh/doc/manuale/node65.html")
        else:
            wavelength = kwargs['wavelength']

        # If wavelength is a single value GenericClient will have made it a
        # Quantity in the kwargs.
        if not isinstance(wavelength, u.Quantity):
            raise ValueError("Wavelength to NORH must be one value not {}.".format(wavelength))

        wavelength = wavelength.to(u.GHz, equivalencies=u.spectral())
        if wavelength == 34 * u.GHz:
            freq = 'tcz'
        elif wavelength == 17 * u.GHz:
            freq = 'tca'
        else:
            raise ValueError("NORH Data can be downloaded for 17GHz or 34GHz,"
                             " see http://solar.nro.nao.ac.jp/norh/doc/manuale/node65.html")

        # If start of time range is before 00:00, converted to such, so
        # files of the requested time ranger are included.
        # This is done because the archive contains daily files.
        if timerange.start.time() != datetime.time(0, 0):
            timerange = TimeRange('{:%Y-%m-%d}'.format(timerange.start),
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
            times.append(TimeRange(t0, t0 + datetime.timedelta(days=1)))
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
        required = {a.Time, a.Instrument}
        optional = {a.Wavelength}
        all_attrs = {type(x) for x in query}

        ops = all_attrs - required
        # If ops is empty or equal to optional we are ok, otherwise we don't
        # match
        if ops and ops != optional:
            return False

        # if we get this far we have either Instrument and Time
        # or Instrument, Time and Wavelength
        for x in query:
            if isinstance(x, a.Instrument) and x.value.lower() == 'norh':
                return True

        return False
