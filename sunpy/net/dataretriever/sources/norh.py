#  Author: Rishabh Sharma <rishabh.sharma.gunner@gmail.com>
#  This Module was developed under funding provided by
#  Google Summer of Code 2014

import astropy.units as u

from sunpy.extern.six.moves.urllib.parse import urljoin

from sunpy.net import attrs as a
from ..client import GenericClient

__all__ = ['NoRHClient']


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
        days = timerange.get_dates()
        urls = []
        for day in days:
            urls.append(self._get_url_for_date(day, **kwargs))
        return urls

    def _get_url_for_date(self, date, **kwargs):
        """
        Return URL for corresponding date.

        Parameters
        ----------
        date : Python datetime object

        Returns
        -------
        string
            The URL for the corresponding date.
        """

        # default urllib password anonymous@ is not accepted by the NoRH FTP
        # server. include an accepted password in base url
        baseurl = 'ftp://anonymous:mozilla@example.com@solar-pub.nao.ac.jp/pub/nsro/norh/data/tcx/'

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
            final_url = urljoin(baseurl,
                                date.strftime('%Y/%m/' + 'tcz' + '%y%m%d'))
        elif wavelength == 17 * u.GHz:
            final_url = urljoin(baseurl,
                                date.strftime('%Y/%m/' + 'tca' + '%y%m%d'))
        else:
            raise ValueError("NORH Data can be downloaded for 17GHz or 34GHz,"
                             " see http://solar.nro.nao.ac.jp/norh/doc/manuale/node65.html")

        return final_url

    def _makeimap(self):
        """
        Helper Function used to hold information about source.
        """
        self.map_['source'] = 'NAOJ'
        self.map_['provider'] = 'NRO'
        self.map_['instrument'] = 'NORH'
        self.map_['phyobs'] = ''

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
