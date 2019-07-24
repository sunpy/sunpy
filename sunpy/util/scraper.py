"""
This module provides a web scraper.
"""
import os
import re
import datetime
from ftplib import FTP
from urllib.error import HTTPError
from urllib.request import urlopen

from bs4 import BeautifulSoup

import astropy.units as u
from astropy.time import Time, TimeDelta

__all__ = ['Scraper']

# regular expressions to convert datetime format
# added `%e` as for milliseconds `%f/1000`
TIME_CONVERSIONS = {'%Y': r'\d{4}', '%y': r'\d{2}',
                    '%b': '[A-Z][a-z]{2}', '%B': r'\W', '%m': r'\d{2}',
                    '%d': r'\d{2}', '%j': r'\d{3}',
                    '%H': r'\d{2}', '%I': r'\d{2}',
                    '%M': r'\d{2}',
                    '%S': r'\d{2}', '%e': r'\d{3}', '%f': r'\d{6}'}


class Scraper:
    """
    A Scraper to scrap web data archives based on dates.

    Parameters
    ----------
    pattern : `str`
        A string containing the url with the date encoded as datetime formats,
        and any other parameter as ``kwargs`` as a string format.

    Attributes
    ----------
    pattern : `str`
        A converted string with the kwargs.
    now : `datetime.datetime`
        The pattern with the actual date.

    Examples
    --------
    >>> # Downloading data from SolarMonitor.org
    >>> from sunpy.util.scraper import Scraper
    >>> solmon_pattern = ('http://solarmonitor.org/data/'
    ...                   '%Y/%m/%d/fits/{instrument}/'
    ...                   '{instrument}_{wave:05d}_fd_%Y%m%d_%H%M%S.fts.gz')
    >>> solmon = Scraper(solmon_pattern, instrument = 'swap', wave = 174)
    >>> print(solmon.pattern)
    http://solarmonitor.org/data/%Y/%m/%d/fits/swap/swap_00174_fd_%Y%m%d_%H%M%S.fts.gz
    >>> print(solmon.now)  # doctest: +SKIP
    http://solarmonitor.org/data/2017/11/20/fits/swap/swap_00174_fd_20171120_193933.fts.gz

    Notes
    -----
    The ``now`` attribute does not return an existent file, but just how the
    pattern looks with the actual time.
    """
    def __init__(self, pattern, **kwargs):
        self.pattern = pattern.format(**kwargs)
        milliseconds = re.search(r'\%e', self.pattern)
        if not milliseconds:
            self.now = datetime.datetime.now().strftime(self.pattern)
        else:
            now = datetime.datetime.now()
            milliseconds_ = int(now.microsecond / 1000.)
            self.now = now.strftime('{start}{milli:03d}{end}'.format(
                start=self.pattern[0:milliseconds.start()],
                milli=milliseconds_,
                end=self.pattern[milliseconds.end():]
            ))

    def matches(self, filepath, date):
        return date.strftime(self.pattern) == filepath

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
            `List` of all the possible directories valid for the time range given.
            Notice that these directories may not exist in the archive.
        """
        # find directory structure - without file names
        directorypattern = os.path.dirname(self.pattern) + '/'
        # TODO what if there's not slashes?
        rangedelta = timerange.dt
        timestep = self._smallerPattern(directorypattern)
        if timestep is None:
            return [directorypattern]
        else:
            # Number of elements in the time range (including end)
            n_steps = rangedelta.sec/timestep.sec
            TotalTimeElements = int(round(n_steps)) + 1
            directories = [(timerange.start + n * timestep).strftime(directorypattern)
                           for n in range(TotalTimeElements)]  # TODO if date <= endate
            return directories

    def _URL_followsPattern(self, url):
        """
        Check whether the url provided follows the pattern.
        """
        pattern = self.pattern
        for k, v in TIME_CONVERSIONS.items():
            pattern = pattern.replace(k, v)
        matches = re.match(pattern, url)
        if matches:
            return matches.end() == matches.endpos == len(self.now)
        return False

    def _extractDateURL(self, url):
        """
        Extracts the date from a particular url following the pattern.
        """
        # remove the user and passwd from files if there:
        url = url.replace("anonymous:data@sunpy.org@", "")

        def url_to_list(txt):
            # Substitutes '.' and '_' for '/'.
            return re.sub(r'\.|_', '/', txt).split('/')

        # create a list of all the blocks in times - assuming they are all
        # separated with either '.', '_' or '/'.
        pattern_list = url_to_list(self.pattern)
        url_list = url_to_list(url)
        time_order = ['%Y', '%y', '%b', '%B', '%m', '%d', '%j',
                      '%H', '%I', '%M', '%S', '%e', '%f']
        final_date = []
        final_pattern = []
        # Find in directory and filename
        for pattern_elem, url_elem in zip(pattern_list, url_list):
            time_formats = [x for x in time_order if x in pattern_elem]
            if len(time_formats) > 0:
                # Find whether there's text that should not be here
                toremove = re.split('%.', pattern_elem)
                if len(toremove) > 0:
                    for bit in toremove:
                        if bit != '':
                            url_elem = url_elem.replace(bit, '', 1)
                            pattern_elem = pattern_elem.replace(bit, '', 1)
                final_date.append(url_elem)
                final_pattern.append(pattern_elem)
                for time_bit in time_formats:
                    time_order.remove(time_bit)
        # Find and remove repeated elements eg: %Y in ['%Y', '%Y%m%d']
        # Make all as single strings
        date_together = ''.join(final_date)
        pattern_together = ''.join(final_pattern)
        re_together = pattern_together
        for k, v in TIME_CONVERSIONS.items():
            re_together = re_together.replace(k, v)

        # Lists to contain the unique elements of the date and the pattern
        final_date = list()
        final_pattern = list()
        re_together = re_together.replace('[A-Z]', '\\[A-Z]')
        for p, r in zip(pattern_together.split('%')[1:], re_together.split('\\')[1:]):
            if p == 'e':
                continue
            regexp = fr'\{r}' if not r.startswith('[') else r
            pattern = f'%{p}'
            date_part = re.search(regexp, date_together)
            date_together = date_together[:date_part.start()] \
                + date_together[date_part.end():]
            if pattern not in final_pattern:
                final_pattern.append(f'%{p}')
                final_date.append(date_part.group())
        return Time.strptime(' '.join(final_date),
                             ' '.join(final_pattern))

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
        >>> from sunpy.util.scraper import Scraper
        >>> solmon_pattern = ('http://solarmonitor.org/data/'
        ...                   '%Y/%m/%d/fits/{instrument}/'
        ...                   '{instrument}_{wave:05d}_fd_%Y%m%d_%H%M%S.fts.gz')
        >>> solmon = Scraper(solmon_pattern, instrument = 'swap', wave = 174)
        >>> from sunpy.time import TimeRange
        >>> timerange = TimeRange('2015-01-01','2015-01-01T16:00:00')
        >>> print(solmon.filelist(timerange))  # doctest: +REMOTE_DATA
        ['http://solarmonitor.org/data/2015/01/01/fits/swap/swap_00174_fd_20150101_025423.fts.gz',
         'http://solarmonitor.org/data/2015/01/01/fits/swap/swap_00174_fd_20150101_061145.fts.gz',
         'http://solarmonitor.org/data/2015/01/01/fits/swap/swap_00174_fd_20150101_093037.fts.gz',
         'http://solarmonitor.org/data/2015/01/01/fits/swap/swap_00174_fd_20150101_124927.fts.gz']

        Note
        ----

        The search is strict with the time range, so if the archive scraped contains daily files,
        but the range doesn't start from the beginning of the day, then the file for that day
        won't be selected. The end of the timerange will normally be OK as includes the file
        on such end time.
        """
        directories = self.range(timerange)
        filesurls = []
        if directories[0][0:3] == "ftp":  # TODO use urlsplit from pr #1807
            return self._ftpfileslist(timerange)
        for directory in directories:
            try:
                opn = urlopen(directory)
                try:
                    soup = BeautifulSoup(opn, "html.parser")
                    for link in soup.find_all("a"):
                        href = link.get("href")
                        if href.endswith(self.pattern.split('.')[-1]):
                            fullpath = directory + href
                            if self._URL_followsPattern(fullpath):
                                datehref = self._extractDateURL(fullpath)
                                if (datehref >= timerange.start and
                                        datehref <= timerange.end):
                                    filesurls.append(fullpath)
                finally:
                    opn.close()
            except HTTPError as http_err:
                # Ignore missing directories (issue #2684).
                if http_err.code == 404:
                    continue
                raise
            except Exception:
                raise
        return filesurls

    def _ftpfileslist(self, timerange):
        directories = self.range(timerange)
        filesurls = list()
        domain = directories[0].find('//')
        domain_slash = directories[0].find('/', 6)  # TODO: Use also urlsplit from pr #1807
        ftpurl = directories[0][domain + 2:domain_slash]
        with FTP(ftpurl, user="anonymous", passwd="data@sunpy.org") as ftp:
            for directory in directories:
                ftp.cwd(directory[domain_slash:])
                for file_i in ftp.nlst():
                    fullpath = directory + file_i
                    if self._URL_followsPattern(fullpath):
                        datehref = self._extractDateURL(fullpath)
                        if (datehref >= timerange.start and
                                datehref <= timerange.end):
                            filesurls.append(fullpath)
        filesurls = ['ftp://anonymous:data@sunpy.org@' + url[domain + 2:] for url in filesurls]
        return filesurls

    def _smallerPattern(self, directoryPattern):
        """
        Obtain the smaller time step for the given pattern.
        """
        try:
            if "%S" in directoryPattern:
                return TimeDelta(1*u.second)
            elif "%M" in directoryPattern:
                return TimeDelta(1*u.minute)
            elif any(hour in directoryPattern for hour in ["%H", "%I"]):
                return TimeDelta(1*u.hour)
            elif any(day in directoryPattern for day in ["%d", "%j"]):
                return TimeDelta(1*u.day)
            elif any(month in directoryPattern for month in ["%b", "%B", "%m"]):
                return TimeDelta(31*u.day)
            elif any(year in directoryPattern for year in ["%Y", "%y"]):
                return TimeDelta(365*u.day)
            else:
                return None
        except Exception:
            raise
