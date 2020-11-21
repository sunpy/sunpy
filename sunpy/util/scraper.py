"""
This module provides a web scraper.
"""
import os
import re
import calendar
import datetime
import warnings
from time import sleep
from ftplib import FTP
from urllib.error import HTTPError
from urllib.parse import urlsplit
from urllib.request import urlopen

from bs4 import BeautifulSoup

import astropy.units as u
from astropy.time import Time, TimeDelta

from sunpy import log
from sunpy.extern.parse import parse
from sunpy.time import TimeRange
from sunpy.util.exceptions import SunpyUserWarning

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

    regex : `bool`
        Set to `True` if parts of the pattern uses regexp symbols. Be careful that
        periods `.` matches any character and therefore it's better to escape them.
        If `regexp` is used, other ``kwargs`` are ignored and string replacement is
        not possible. Default is `False`.

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
    def __init__(self, pattern, regex=False, **kwargs):
        if regex:
            self.pattern = pattern
            if kwargs:
                warnings.warn('regexp being used, the extra arguments passed are being ignored',
                              SunpyUserWarning)
        else:
            self.pattern = pattern.format(**kwargs)
        self.domain = "{0.scheme}://{0.netloc}/".format(urlsplit(self.pattern))
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
            return matches.end() == matches.endpos
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
        if urlsplit(directories[0]).scheme == "ftp":
            return self._ftpfileslist(timerange)
        if urlsplit(directories[0]).scheme == "file":
            return self._localfilelist(timerange)
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
                            if self._URL_followsPattern(fullpath):
                                if self._check_timerange(fullpath, timerange):
                                    filesurls.append(fullpath)
                finally:
                    opn.close()
            except HTTPError as http_err:
                # Ignore missing directories (issue #2684).
                if http_err.code == 404:
                    continue
                if http_err.code == 429:
                    # See if the server has told us how long to back off for
                    retry_after = http_err.hdrs.get('Retry-After', 1)
                    try:
                        # Ensure that we can parse the header as an int in sec
                        retry_after = int(retry_after)
                    except Exception:
                        retry_after = 1

                    log.debug(f"Got 429 while scraping {directory}, waiting for {retry_after} seconds before retrying.")

                    sleep(retry_after)

                    # Put this dir back on the queue
                    directories.insert(0, directory)
                    continue
                raise
            except Exception:
                raise
        return filesurls

    def _ftpfileslist(self, timerange):
        directories = self.range(timerange)
        filesurls = list()
        ftpurl = urlsplit(directories[0]).netloc
        with FTP(ftpurl, user="anonymous", passwd="data@sunpy.org") as ftp:
            for directory in directories:
                ftp.cwd(urlsplit(directory).path)
                for file_i in ftp.nlst():
                    fullpath = directory + file_i
                    if self._URL_followsPattern(fullpath):
                        if self._check_timerange(fullpath, timerange):
                            filesurls.append(fullpath)

        filesurls = [f'ftp://' + "{0.netloc}{0.path}".format(urlsplit(url))
                     for url in filesurls]

        return filesurls

    def _localfilelist(self, timerange):
        pattern = self.pattern
        pattern_temp = pattern.replace('file://', '')
        if os.name == 'nt':
            pattern_temp = pattern_temp.replace('\\', '/')
            prefix = 'file:///'
        else:
            prefix = 'file://'
        self.pattern = pattern_temp
        directories = self.range(timerange)
        filepaths = list()
        for directory in directories:
            for file_i in os.listdir(directory):
                fullpath = directory + file_i
                if self._URL_followsPattern(fullpath):
                    if self._check_timerange(fullpath, timerange):
                        filepaths.append(fullpath)
        filepaths = [prefix + path for path in filepaths]
        self.pattern = pattern
        return filepaths

    def _check_timerange(self, url, timerange):
        """
        Checks whether the time extracted from the URL
        is valid according to the given time range.

        Parameters
        ----------
        url: `str`
            URL of the file.
        timerange : `~sunpy.time.TimeRange`
            Time interval for which files were searched.

        Returns
        -------
        `bool`
            `True` if URL's time overlaps the given timerange, else `False`.
        """
        if hasattr(self, 'extractor'):
            exdict = parse(self.extractor, url).named
            tr = get_timerange_from_exdict(exdict)
            return (tr.end >= timerange.start and tr.start <= timerange.end)
        else:
            datehref = self._extractDateURL(url).to_datetime()
            return (timerange.start.to_datetime() <= datehref <= timerange.end.to_datetime())

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

    def _extract_files_meta(self, timerange, extractor, matcher=None):
        """
        Returns metadata information contained in URLs.

        Parameters
        ----------
        timerange : `~sunpy.time.TimeRange`
            Time interval where to find the directories for a given pattern.
        extractor: `str`
            Pattern to extract metadata by parsing the URL.
        matcher: `dict`
            Dictionary to check if extracted metadata is valid.

        Returns
        -------
        `list` of `dict`
            List of metadata info for all URLs.
        """
        self.extractor = extractor
        urls = self.filelist(timerange)
        metalist = []
        for url in urls:
            metadict = parse(extractor, url)
            if metadict is not None:
                append = True
                metadict = metadict.named
                metadict['url'] = url
                if matcher is not None:
                    for k in metadict:
                        if k in matcher and str(metadict[k]) not in matcher[k]:
                            append = False
                            break
                if append:
                    metalist.append(metadict)
        return metalist


def get_timerange_from_exdict(exdict):
    """
    Function to get URL's timerange using extracted metadata.
    It computes start and end times first using the given
    dictionary and then returns a timerange.

    Parameters
    ----------
    exdict : `dict`
        Metadata extracted from the file's url.

    Returns
    -------
    `~sunpy.time.TimeRange`
        The time range of the file.
    """
    datetypes = ['year', 'month', 'day']
    timetypes = ['hour', 'minute', 'second', 'millisecond']
    dtlist = [int(exdict.get(d, 1)) for d in datetypes]
    dtlist.extend([int(exdict.get(t, 0)) for t in timetypes])
    startTime = Time(datetime.datetime(*dtlist))

    tdelta = 1*u.millisecond
    if "year" in exdict:
        if calendar.isleap(int(exdict['year'])):
            tdelta = 366*u.day
        else:
            tdelta = 365*u.day
    if "month" in exdict:
        days_in_month = calendar.monthrange(int(exdict['year']), int(exdict['month']))[1]
        tdelta = days_in_month*u.day
    if "day" in exdict:
        tdelta = 1*u.day
    if "hour" in exdict:
        tdelta = 1*u.hour
    if "minute" in exdict:
        tdelta = 1*u.minute
    if "second" in exdict:
        tdelta = 1*u.second

    endTime = startTime + TimeDelta(tdelta) - TimeDelta(1*u.millisecond)
    file_timerange = TimeRange(startTime, endTime)
    return file_timerange
