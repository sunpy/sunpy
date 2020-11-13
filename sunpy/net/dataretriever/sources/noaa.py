# Author: Rishabh Sharma <rishabh.sharma.gunner@gmail.com>
# This module was developed under funding provided by
# Google Summer of Code 2014
import pathlib
import tarfile
from collections import OrderedDict

from astropy import units as u
from astropy.time import Time, TimeDelta

from sunpy.net.dataretriever import GenericClient
from sunpy.util.parfive_helpers import Downloader

__all__ = ['NOAAIndicesClient', 'NOAAPredictClient', 'SRSClient']


class NOAAIndicesClient(GenericClient):
    """
    Provides access to the NOAA solar cycle indices.

    Uses the `SWPC NOAA archive <https://services.swpc.noaa.gov/json/solar-cycle/>`__.
    This is a fixed dataset so the result is independent of the time range.

    Examples
    --------

    >>> from sunpy.net import Fido, attrs as a
    >>> results = Fido.search(a.Time("2016/1/1", "2016/1/2"),
    ...                       a.Instrument.noaa_indices)  #doctest: +REMOTE_DATA
    >>> results  #doctest: +REMOTE_DATA
    <sunpy.net.fido_factory.UnifiedResponse object at ...>
    Results from 1 Provider:
    <BLANKLINE>
    1 Results from the NOAAIndicesClient:
         Start Time           End Time      Source  Instrument  Wavelength
    ------------------- ------------------- ------ ------------ ----------
    2016-01-01 00:00:00 2016-01-02 00:00:00   sdic noaa-indices        nan
    <BLANKLINE>
    <BLANKLINE>

    """
    @staticmethod
    def _get_url_for_timerange(timerange, **kwargs):
        """
        Helper function:
        """
        return ["https://services.swpc.noaa.gov/json/solar-cycle/observed-solar-cycle-indices.json"]

    def _makeimap(self):
        """
        Helper Function:used to hold information about source.
        """
        self.map_['source'] = 'sdic'
        self.map_['instrument'] = 'noaa-indices'
        self.map_['physobs'] = 'sunspot number'
        self.map_['provider'] = 'swpc'

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
        chkattr = ['Time', 'Instrument']
        chklist = [x.__class__.__name__ in chkattr for x in query]
        for x in query:
            if x.__class__.__name__ == 'Instrument' and x.value == 'noaa-indices':
                return all(chklist)
        return False

    @classmethod
    def register_values(cls):
        from sunpy.net import attrs
        adict = {attrs.Instrument: [
            ('NOAA-Indices', 'Recent Solar Indices of Observed Monthly Mean Values')]}
        return adict


class NOAAPredictClient(GenericClient):
    """
    Provides access to the NOAA SWPC predicted sunspot Number and 10.7 cm radio flux values.

    Uses the `SWPC NOAA archive <https://services.swpc.noaa.gov/json/solar-cycle/>`__.
    This is a fixed prediction so the result is independent of the time range.

    Examples
    --------

    >>> from sunpy.net import Fido, attrs as a
    >>> results = Fido.search(a.Time("2016/1/1", "2016/1/2"),
    ...                       a.Instrument.noaa_predict)  #doctest: +REMOTE_DATA
    >>> results  #doctest: +REMOTE_DATA
    <sunpy.net.fido_factory.UnifiedResponse object at ...>
    Results from 1 Provider:
    <BLANKLINE>
    1 Results from the NOAAPredictClient:
         Start Time           End Time      Source  Instrument  Wavelength
    ------------------- ------------------- ------ ------------ ----------
    2016-01-01 00:00:00 2016-01-02 00:00:00   ises noaa-predict        nan
    <BLANKLINE>
    <BLANKLINE>

    """
    @staticmethod
    def _get_default_uri():
        """Return the url to download indices"""
        return ["https://services.swpc.noaa.gov/json/solar-cycle/predicted-solar-cycle.json"]

    def _get_url_for_timerange(self, timerange, **kwargs):
        """
        Helper function:
        """
        return NOAAPredictClient._get_default_uri()

    def _makeimap(self):
        """
        Helper Function:used to hold information about source.
        """
        self.map_['source'] = 'ises'
        self.map_['instrument'] = 'noaa-predict'
        self.map_['physobs'] = 'sunspot number'
        self.map_['provider'] = 'swpc'

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
        chkattr = ['Time', 'Instrument']
        chklist = [x.__class__.__name__ in chkattr for x in query]
        for x in query:
            if x.__class__.__name__ == 'Instrument' and x.value.lower() == 'noaa-predict':
                return all(chklist)
        return False

    @classmethod
    def register_values(cls):
        from sunpy.net import attrs
        adict = {attrs.Instrument: [
            ('NOAA-Predict', 'Predicted Sunspot Number And Radio Flux Values With Expected Ranges.')]}
        return adict


class SRSClient(GenericClient):
    """
    Provides access to the NOAA SWPC solar region summary data.

    Uses the `ftp archive <ftp://ftp.swpc.noaa.gov/pub/warehouse/>`__.

    Notes
    -----
    Data pre-1996 is in free-form text, which cannot be parsed by sunpy, and
    therefore only results from 1996 onwards are returned by this client.

    Examples
    --------
    >>> from sunpy.net import Fido, attrs as a
    >>> results = Fido.search(a.Time("2016/1/1", "2016/1/2"),
    ...                       a.Instrument.soon)  #doctest: +REMOTE_DATA
    >>> results  #doctest: +REMOTE_DATA
    <sunpy.net.fido_factory.UnifiedResponse object at ...>
    Results from 1 Provider:
    <BLANKLINE>
    2 Results from the SRSClient:
         Start Time           End Time        Source  Instrument Wavelength
    ------------------- ------------------- --------- ---------- ----------
    2016-01-01 00:00:00 2016-01-02 00:00:00 NOAA/USAF       SOON        nan
    2016-01-01 00:00:00 2016-01-02 00:00:00 NOAA/USAF       SOON        nan
    <BLANKLINE>
    <BLANKLINE>

    """

    def _get_url_for_timerange(self, timerange, **kwargs):
        result = list()
        base_url = 'ftp://ftp.swpc.noaa.gov/pub/warehouse/'
        total_days = int(timerange.days.value) + 1
        all_dates = timerange.split(total_days)
        today_year = int(Time.now().strftime('%Y'))
        for day in all_dates:
            end_year = int(day.end.strftime('%Y'))
            if end_year > today_year or end_year < 1996:
                continue
            elif end_year == today_year:
                suffix = '{}/SRS/{}SRS.txt'.format(
                    end_year, day.end.strftime('%Y%m%d'))
            else:
                suffix = '{}/{}_SRS.tar.gz'.format(
                    end_year, day.end.strftime('%Y'))
            url = base_url + suffix
            result.append(url)
        return result

    def fetch(self, qres, path=None, error_callback=None, **kwargs):
        """
        Download a set of results.

        Parameters
        ----------
        qres : `~sunpy.net.dataretriever.QueryResponse`
            Results to download.

        Returns
        -------
        Results Object
        """

        urls = [qrblock.url for qrblock in qres]

        filenames = []
        local_filenames = []

        for i, [url, qre] in enumerate(zip(urls, qres)):
            name = url.split('/')[-1]

            # temporary fix !!! coz All QRBs have same start_time values
            day = Time(qre.time.start.strftime('%Y-%m-%d')) + TimeDelta(i*u.day)

            if name not in filenames:
                filenames.append(name)

            if name.endswith('.gz'):
                local_filenames.append('{}SRS.txt'.format(day.strftime('%Y%m%d')))
            else:
                local_filenames.append(name)

        if path is not None:
            path = pathlib.Path(path)
        # Files to be actually downloaded
        paths = self._get_full_filenames(qres, filenames, path)

        # Those files that will be present after get returns
        local_paths = self._get_full_filenames(qres, local_filenames, path)

        # remove duplicate urls. This will make paths and urls to have same number of elements.
        # OrderedDict is required to maintain ordering because it will be zipped with paths later
        urls = list(OrderedDict.fromkeys(urls))

        dobj = Downloader(max_conn=5)

        for aurl, fname in zip(urls, paths):
            dobj.enqueue_file(aurl, filename=fname)

        paths = dobj.download()

        outfiles = []
        for fname, srs_filename in zip(local_paths, local_filenames):

            name = fname.name

            past_year = False
            for i, fname2 in enumerate(paths):
                fname2 = pathlib.Path(fname2)

                if fname2.name.endswith('.txt'):
                    continue

                year = fname2.name.split('_SRS')[0]

                if year in name:
                    TarFile = tarfile.open(fname2)
                    filepath = fname.parent
                    member = TarFile.getmember('SRS/' + srs_filename)
                    member.name = name
                    TarFile.extract(member, path=filepath)
                    TarFile.close()

                    outfiles.append(fname)

                    past_year = True
                    break

            if past_year is False:
                outfiles.append(fname)

        paths.data = list(map(str, outfiles))
        return paths

    def _makeimap(self):
        self.map_['source'] = 'swpc'
        self.map_['instrument'] = 'SOON'
        self.map_['physobs'] = 'SRS'
        self.map_['source'] = 'NOAA/USAF'

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
        from sunpy.net import attrs as a
        if a.Time not in [type(at) for at in query]:
            return False
        for x in query:
            if (x.__class__.__name__ == "Instrument" and
                    str(x.value).lower() in ["soon", "srs-table"]):
                return True
        return False

    @classmethod
    def register_values(cls):
        from sunpy.net import attrs
        adict = {attrs.Instrument: [("SOON", "Solar Region Summary."),
                                    ("SRS-Table", "Solar Region Summary.")]}
        return adict
