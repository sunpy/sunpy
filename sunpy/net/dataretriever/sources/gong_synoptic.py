import re
import os
import gzip
from pathlib import Path

from urllib.error import HTTPError
from urllib.request import urlopen

from bs4 import BeautifulSoup

from parfive import Downloader

from ..client import GenericClient
from sunpy.util.scraper import Scraper
from sunpy.time import parse_time

__all__ = ['GongSynopticClient']


class GongSynopticClient(GenericClient):
    """
    Provides access to the merged data product for NSO-GONG synoptic maps
    `archive <gong2.nso.edu/oQR/zqs/>`__

    Examples
    --------

    >>> from sunpy.net import Fido, attrs as a
    >>> results = Fido.search(a.Time("2019/12/31 22:00", "2020/1/1 02:00"),
    ...                       a.Instrument('GONG'))  #doctest: +REMOTE_DATA
    >>> results  #doctest: +REMOTE_DATA
    <sunpy.net.fido_factory.UnifiedResponse object at ...>
    Results from 1 Provider:
    <BLANKLINE>
    3 Results from the GongSynopticClient:
         Start Time           End Time      Source Instrument Wavelength
    ------------------- ------------------- ------ ---------- ----------
    2019-12-31 22:00:00 2020-01-01 02:00:00    NSO       GONG        nan
    2019-12-31 22:00:00 2020-01-01 02:00:00    NSO       GONG        nan
    2019-12-31 22:00:00 2020-01-01 02:00:00    NSO       GONG        nan
    <BLANKLINE>
    <BLANKLINE>

    """

    def _get_url_for_timerange(self, timerange, **kwargs):
        """
        Return URL(s) for corresponding timerange.

        Parameters
        ----------
        timerange : `~sunpy.time.TimeRange`
            The time range you want the files for.

        Returns
        -------
        `list`
            The URL(s) for the corresponding timerange.
        """
        prefix = 'https://gong2.nso.edu/oQR/zqs/%Y%m/mrzqs%y%m%d/'
        file_pattern = 'mrzqs\d{6}t\d{4}c\d+_\d+.fits.gz'
        pref_scrap = Scraper(prefix)
        directories = pref_scrap.range(timerange)
        filesurls = []
        for directory in directories:
            try:
                opn = urlopen(directory)
                try:
                    soup = BeautifulSoup(opn, "html.parser")
                    for link in soup.find_all("a"):
                        href = link.get("href")
                        if re.match(file_pattern, href):
                            # extracting timestamp from the filename
                            timestring = href[5:16]
                            tstamp = parse_time('20'+timestring.upper()+'00')
                            if tstamp >= timerange.start and tstamp <= timerange.end:
                                fullpath = directory + href
                                filesurls.append(fullpath)
                finally:
                    opn.close()
            except HTTPError as http_err:
                if http_err.code == 404:
                    continue
                raise
            except Exception:
                raise
        return filesurls

    def fetch(self, qres, path=None, error_callback=None, **kwargs):
        """
        Download a set of results.

        Parameters
        ----------
        qres : `~sunpy.net.dataretriever.QueryResponse`
            Results to download.

        path : `str`, optional
            The directory to retrieve the files into.

        Returns
        -------
        results: `parfive.Results`

        """
        if path is not None:
            path = Path(path)

        urls = [qrblock.url for qrblock in qres]

        filenames = [url.split('/')[-1] for url in urls]
        local_filenames = [name[:-3] for name in filenames]
        paths = self._get_full_filenames(qres, filenames, path)
        local_paths = self._get_full_filenames(qres, local_filenames, path)

        downloader = Downloader()

        for url, filename in zip(urls, paths):
            downloader.enqueue_file(url, filename=filename)

        resobj = downloader.download()

        for loc_fpath, fpath in zip(local_paths, paths):
            f = gzip.open(fpath, 'r')
            g = gzip.open(loc_fpath, 'w')
            g.write(f.read())

        resobj.data = list(map(str, local_paths))
        return resobj

    def _makeimap(self):
        """
        Helper Function:used to hold information about source.
        """
        self.map_['instrument'] = 'GONG'
        self.map_['source'] = 'NSO'

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
            if x.__class__.__name__ == 'Instrument' and x.value.lower() == 'gong':
                return all(chklist)
        return False
