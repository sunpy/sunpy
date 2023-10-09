"""
This module provides a web scraper.
"""
import os
import re
from time import sleep
from ftplib import FTP
from datetime import datetime
from urllib.error import HTTPError
from urllib.parse import urlsplit
from urllib.request import urlopen

from bs4 import BeautifulSoup

from sunpy import log
from sunpy.extern.parse import parse
from sunpy.net.scraper_utils import check_timerange, date_floor, extract_timestep

__all__ = ['Scraper']

# `parse` expressions to convert into datetime format
TIME_CONVERSIONS = {"{year:4d}": "%Y", "{year:2d}": "%y",
            "{month:2d}": "%m",
            "{month_name:l}": "%B",
            "{month_name_abbr:l}": "%b",
            "{day:2d}": "%d", "{day_of_year:3d}": "%j",
            "{hour:2d}": "%H",
            "{minute:2d}": "%M",
            "{second:2d}": "%S",
            "{microsecond:6d}": "%f",
            "{millisecond:3d}": "%e", # added `%e` as for milliseconds `%f/1000`
            "{week_number:2d}": "%W",
        }


class Scraper:
    """
    A Scraper to scrap web data archives based on dates.

    Parameters
    ----------
    pattern : `str`
        A string containing the url with the date and other information to be
        extracted encoded as ``parse`` formats, and any other ``kwargs`` parameters
        as a string format, the former represented using double curly-brackets
        to differentiate from the latter.
        The accepted parse representations for datetime values are as given in ``TIME_CONVERSIONS``.
        This can also be a uri to a local file patterns.


    Attributes
    ----------
    pattern : `str`
        A converted string with the kwargs.
    now : `datetime.datetime`
        The pattern with the actual date.

    Examples
    --------
    >>> from sunpy.net import Scraper
    >>> pattern = ('http://proba2.oma.be/{instrument}/data/bsd/{{year:4d}}/{{month:2d}}/{{day:2d}}/'
    ...            '{instrument}_lv1_{{year:4d}}{{month:2d}}{{day:2d}}_{{hour:2d}}{{month:2d}}{{second:2d}}.fits')
    >>> swap = Scraper(pattern, instrument='swap')
    >>> print(swap.pattern)
    http://proba2.oma.be/swap/data/bsd/{year:4d}/{month:2d}/{day:2d}/swap_lv1_{year:4d}{month:2d}{day:2d}_{hour:2d}{month:2d}{second:2d}.fits
    >>> print(swap.now)  # doctest: +SKIP
    http://proba2.oma.be/swap/data/bsd/2022/12/21/swap_lv1_20221221_112433.fits

    Notes
    -----
    The ``now`` attribute does not return an existent file, but just how the
    pattern looks with the actual time.
    """

    def __init__(self, pattern, **kwargs):
        pattern = pattern.format(**kwargs)
        timepattern = pattern
        for k, v in TIME_CONVERSIONS.items():
            if k in timepattern:
                timepattern = timepattern.replace(k,v)
        self.timepattern = timepattern
        if "year:4d" in pattern and "year:2d" in pattern:
            pattern = pattern.replace("year:2d", ":2d")
        self.pattern = pattern
        self.domain = f"{urlsplit(self.pattern).scheme}://{urlsplit(self.pattern).netloc}/"
        milliseconds = re.search(r'\%e', self.timepattern)
        if not milliseconds:
            self.now = datetime.now().strftime(self.timepattern)
        else:
            now = datetime.now()
            milliseconds_ = int(now.microsecond / 1000.)
            self.now = now.strftime('{start}{milli:03d}{end}'.format(
                start=self.timepattern[0:milliseconds.start()],
                milli=milliseconds_,
                end=self.timepattern[milliseconds.end():]
            ))

    def matches(self, filepath, date):
        """
        Checks if the given filepath is how the file path is expected
        to look on given date based on the pattern.

        Parameters
        ----------
        filepath : `str`
            File path to check.
        date : `datetime.datetime` or `astropy.time.Time`
            The date for which to check.

        Returns
        -------
        `bool`
            `True` if the given filepath matches with the calculated one for given date, else `False`.
        """
        return parse(date.strftime(self.timepattern), filepath) is not None

    def range(self, timerange):
        """
        Gets the directories for a certain range of time.

        Parameters
        ----------
        timerange : `~sunpy.time.timerange.TimeRange`
            Time interval where to find the directories for a given pattern.

        Returns
        -------
        `list` of `str`
            All the possible directories valid for the time range given.
            Notice that these directories may not exist in the archive.
        """
        # find directory structure - without file names
        if '/' in self.timepattern:
            directorypattern = '/'.join(self.timepattern.split('/')[:-1]) + '/'
        timestep = extract_timestep(directorypattern)
        if timestep is None:
            return [directorypattern]
        else:
            directories = []
            cur = date_floor(timerange.start, timestep)
            end = date_floor(timerange.end, timestep) + timestep
            while cur < end:
                directories.append(cur.strftime(directorypattern))
                cur = cur + timestep
            return directories

    def filelist(self, timerange):
        """
        Returns the list of existent files in the archive for the given time
        range.

        Parameters
        ----------
        timerange : `~sunpy.time.TimeRange`
            Time interval where to find the directories for a given pattern.

        Returns
        -------
        filesurls : `list` of `str`
            List of all the files found between the time range given.

        Examples
        --------
        >>> from sunpy.net import Scraper
        >>> pattern = ('http://proba2.oma.be/{instrument}/data/bsd/{{year:4d}}/{{month:2d}}/{{day:2d}}/'
        ...            '{instrument}_lv1_{{year:4d}}{{month:2d}}{{day:2d}}_{{hour:2d}}{{minute:2d}}{{second:2d}}.fits')
        >>> swap = Scraper(pattern, instrument='swap')
        >>> from sunpy.time import TimeRange
        >>> timerange = TimeRange('2015-01-01T00:08:00','2015-01-01T00:12:00')
        >>> print(swap.filelist(timerange))  # doctest: +REMOTE_DATA
        ['http://proba2.oma.be/swap/data/bsd/2015/01/01/swap_lv1_20150101_000857.fits',
         'http://proba2.oma.be/swap/data/bsd/2015/01/01/swap_lv1_20150101_001027.fits',
         'http://proba2.oma.be/swap/data/bsd/2015/01/01/swap_lv1_20150101_001157.fits']

        While writing the pattern, we can also leverage parse capabilities by using the ``{{}}`` notation to match parts of the filename that cannot be known beforehand:
        >>> from sunpy.net import Scraper
        >>> from sunpy.time import TimeRange
        >>> pattern = 'http://proba2.oma.be/lyra/data/bsd/{{year:4d}}/{{month:2d}}/{{day:2d}}/{{}}_lev{{Level:1d}}_std.fits'
        >>> lyra = Scraper(pattern)
        >>> timerange = TimeRange('2023-03-06T00:08:00','2023-03-12T00:12:00')
        >>> print(swap.filelist(timerange)) # doctest: +REMOTE_DATA
        ['http://proba2.oma.be/lyra/data/bsd/2023/03/06/lyra_20230306-000000_lev2_std.fits',
        'http://proba2.oma.be/lyra/data/bsd/2023/03/06/lyra_20230306-000000_lev3_std.fits',
        '...',
        'http://proba2.oma.be/lyra/data/bsd/2023/03/12/lyra_20230312-000000_lev3_std.fits']

        Notes
        -----
        The search is strict with the time range, so if the archive scraped contains daily files,
        but the range doesn't start from the beginning of the day, then the file for that day
        won't be selected. The end of the timerange will normally be OK as includes the file
        on such end time.
        """
        directories = self.range(timerange)
        if urlsplit(directories[0]).scheme == "ftp":
            return self._ftpfilelist(timerange)
        elif urlsplit(directories[0]).scheme == "file":
            return self._localfilelist(timerange)
        elif urlsplit(directories[0]).scheme == "http":
            return self._httpfilelist(timerange)
        else:
            return ValueError("The provided pattern should either be an FTP or a local file-path, or an HTTP address.")

    def _ftpfilelist(self, timerange):
        """
        Goes over archives available over ftp to return list of files in the given timerange.
        """
        directories = self.range(timerange)
        filesurls = list()
        ftpurl = urlsplit(directories[0]).netloc
        with FTP(ftpurl, user="anonymous", passwd="data@sunpy.org") as ftp:
            for directory in directories:
                try:
                    ftp.cwd(urlsplit(directory).path)
                except Exception as e:
                    log.debug(f"FTP CWD: {e}")
                    continue
                for file_i in ftp.nlst():
                    fullpath = directory + file_i
                    if parse(self.pattern, fullpath):
                        if check_timerange(self.pattern, fullpath, timerange):
                            filesurls.append(fullpath)

        filesurls = ['ftp://' + f"{urlsplit(url).netloc}{urlsplit(url).path}"
                     for url in filesurls]

        return filesurls

    def _localfilelist(self, timerange):
        """
        Goes over locally stored archives to return list of files in the given timerange.
        """
        pattern, timepattern = self.pattern, self.timepattern
        pattern_temp, timepattern_temp = pattern.replace('file://', ''), timepattern.replace('file://', '')
        if os.name == 'nt':
            pattern_temp = pattern_temp.replace('\\', '/')
            prefix = 'file:///'
        else:
            prefix = 'file://'
        # Change pattern variables class-wide
        self.pattern, self.timepattern = pattern_temp, timepattern_temp
        directories = self.range(timerange)
        filepaths = list()
        for directory in directories:
            for file_i in os.listdir(directory):
                fullpath = directory + file_i
                if parse(self.pattern, fullpath):
                    if check_timerange(self.pattern, fullpath, timerange):
                        filepaths.append(fullpath)
        filepaths = [prefix + path for path in filepaths]
        # Set them back to their original values
        self.pattern, self.timepattern = pattern, timepattern
        return filepaths

    def _httpfilelist(self, timerange):
        """
        Goes over http archives hosted on the web, to return list of files in the given timerange.
        """
        directories = self.range(timerange)
        filesurls = list()
        while directories:
            directory = directories.pop(0)
            try:
                opn = urlopen(directory)
                try:
                    soup = BeautifulSoup(opn, "html.parser")
                    for link in soup.find_all("a"):
                        href = link.get("href")
                        if href is not None and href.endswith(self.pattern.split('.')[-1]):
                            if href[0] == '/':
                                fullpath = self.domain + href[1:]
                            else:
                                fullpath = directory + href
                            if parse(self.pattern, fullpath):
                                if check_timerange(self.pattern, fullpath, timerange):
                                    filesurls.append(fullpath)
                finally:
                    opn.close()
            except HTTPError as http_err:
                # Ignore missing directories (issue #2684).
                if http_err.code == 404:
                    continue
                if http_err.code == 429:
                    # See if the server has told us how long to back off for
                    retry_after = http_err.hdrs.get('Retry-After', 2)
                    try:
                        # Ensure that we can parse the header as an int in sec
                        retry_after = int(retry_after)
                    except Exception as e:
                        log.debug(f"Converting retry_after failed: {e}")
                        retry_after = 2
                    log.debug(
                        f"Got 429 while scraping {directory}, waiting for {retry_after} seconds before retrying."
                    )
                    sleep(retry_after)
                    # Put this dir back on the queue
                    directories.insert(0, directory)
                    continue
                raise
            except Exception:
                raise
        return filesurls

    def _extract_files_meta(self, timerange, matcher=None):
        """
        Returns metadata information contained in URLs.

        Parameters
        ----------
        timerange : `~sunpy.time.TimeRange`
            Time interval where to find the directories for a given pattern.
        matcher : `dict`
            Dictionary to check if extracted metadata is valid.

        Returns
        -------
        `list` of `dict`
            List of metadata info for all URLs.
        """
        urls = self.filelist(timerange)
        metalist = []
        for url in urls:
            metadict = parse(self.pattern, url)
            if metadict is not None:
                append = True
                metadict = metadict.named
                metadict['url'] = url
                if 'month' not in metadict:
                    if 'month_name' in metadict:
                        metadict['month'] = datetime.strptime(metadict['month_name'], '%B').month
                    elif 'month_name_abbr' in metadict:
                        metadict['month'] = datetime.strptime(metadict['month_name_abbr'], '%b').month
                if matcher is not None:
                    for k in metadict:
                        if k in matcher and str(metadict[k]) not in matcher[k]:
                            append = False
                            break
                if append:
                    metalist.append(metadict)
        return metalist
