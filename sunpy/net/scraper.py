"""
This module provides a web scraper.
"""
import os
import re
import calendar
from time import sleep
from ftplib import FTP
from datetime import datetime, timedelta
from urllib.error import HTTPError
from urllib.parse import urlsplit
from urllib.request import urlopen

from bs4 import BeautifulSoup
from dateutil.relativedelta import relativedelta

from sunpy import log
from sunpy.extern.parse import parse
from sunpy.time import TimeRange

__all__ = ['Scraper']

# parse expressions to convert into datetime format
# added `%e` as for milliseconds `%f/1000`
TIME_CONVERSIONS = {"{year:4d}": "%Y", "{year:2d}": "%y",
            "{month:2d}": "%m",
            "{month_name:l}": "%B",
            "{month_name_abbr:l}": "%b",
            "{day:2d}": "%d", "{day_of_year:3d}": "%j",
            "{hour:2d}": "%H",
            "{minute:2d}": "%M",
            "{second:2d}": "%S",
            "{microsecond:6d}": "%f",
            "{millisecond:3d}": "%e",
            "{week_number:2d}": "%W",
        }
TIME_QUANTITIES = {'day': timedelta(days=1),
                   'hour': timedelta(hours=1),
                   'minute': timedelta(minutes=1),
                   'second': timedelta(seconds=1),
                   'millisecond': timedelta(milliseconds=1)}


class Scraper:
    """
    A Scraper to scrap web data archives based on dates.

    Parameters
    ----------
    pattern : `str`
        A string containing the url with the date encoded as datetime formats,
        and any other parameter as ``kwargs`` as a string format.
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
        self.domain = "{0.scheme}://{0.netloc}/".format(urlsplit(self.pattern))
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
        return parse(date.strftime(self.timepattern), filepath)

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
        timestep = self._smallerPattern(directorypattern)
        if timestep is None:
            return [directorypattern]
        else:
            directories = []
            cur = self._date_floor(timerange.start, timestep)
            end = self._date_floor(timerange.end, timestep) + timestep
            while cur < end:
                directories.append(cur.strftime(directorypattern))
                cur = cur + timestep
            return directories

    @staticmethod
    def _date_floor(date, timestep):
        """
        Return the "floor" of the given date and time step.

        Parameters
        ----------
        datetime : `datetime.datetime` or `astropy.time.Time`
            The date to floor
        timestep : `dateutil.relativedelta.relativedelta`
            The smallest time step to floor
        Returns
        -------
        `datetime.datetime`
            The time floored at the given time step

        """
        date_parts = [int(p) for p in date.strftime('%Y,%m,%d,%H,%M,%S').split(',')]
        date_parts[-1] = date_parts[-1] % 60
        date = datetime(*date_parts)
        orig_time_tup = date.timetuple()
        time_tup = [orig_time_tup.tm_year, orig_time_tup.tm_mon, orig_time_tup.tm_mday,
                    orig_time_tup.tm_hour, orig_time_tup.tm_min, orig_time_tup.tm_sec]
        if timestep == relativedelta(minutes=1):
            time_tup[-1] = 0
        elif timestep == relativedelta(hours=1):
            time_tup[-2:] = [0, 0]
        elif timestep == relativedelta(days=1):
            time_tup[-3:] = [0, 0, 0]
        elif timestep == relativedelta(months=1):
            time_tup[-4:] = [1, 0, 0, 0]
        elif timestep == relativedelta(years=1):
            time_tup[-5:] = [1, 1, 0, 0, 0]

        return datetime(*time_tup)

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

        Notes
        -----
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
                            if parse(self.pattern, fullpath):
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

    def _ftpfileslist(self, timerange):
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
                        if self._check_timerange(fullpath, timerange):
                            filesurls.append(fullpath)

        filesurls = ['ftp://' + "{0.netloc}{0.path}".format(urlsplit(url))
                     for url in filesurls]

        return filesurls

    def _localfilelist(self, timerange):
        pattern = self.pattern
        pattern_temp = pattern.replace('file://', '')
        timepattern = self.timepattern
        timepattern_temp = timepattern.replace('file://', '')
        if os.name == 'nt':
            pattern_temp = pattern_temp.replace('\\', '/')
            prefix = 'file:///'
        else:
            prefix = 'file://'
        self.pattern = pattern_temp
        self.timepattern = timepattern_temp
        directories = self.range(timerange)
        filepaths = list()
        for directory in directories:
            for file_i in os.listdir(directory):
                fullpath = directory + file_i
                if parse(self.pattern, fullpath):
                    if self._check_timerange(fullpath, timerange):
                        filepaths.append(fullpath)
        filepaths = [prefix + path for path in filepaths]
        self.pattern = pattern
        self.timepattern = timepattern
        return filepaths

    def _check_timerange(self, url, timerange):
        """
        Checks whether the time range represented in *url* intersects
        with the given time range.

        Parameters
        ----------
        url : `str`
            URL of the file.
        timerange : `~sunpy.time.TimeRange`
            Time interval for which files were searched.

        Returns
        -------
        `bool`
            `True` if URL's valid time range overlaps the given timerange, else `False`.
        """
        exdict = parse(self.pattern, url).named
        if exdict['year'] < 100:
            exdict['year'] = 2000 + exdict['year']
        if 'month' not in exdict:
                    if 'month_name' in exdict:
                        exdict['month'] = datetime.strptime(exdict['month_name'], '%B').month
                    elif 'month_name_abbr' in exdict:
                        exdict['month'] = datetime.strptime(exdict['month_name_abbr'], '%b').month
        tr = get_timerange_from_exdict(exdict)
        return tr.intersects(timerange)

    def _smallerPattern(self, directoryPattern):
        """
        Obtain the smaller time step for the given pattern.
        """
        try:
            if "%S" in directoryPattern:
                return relativedelta(seconds=1)
            elif "%M" in directoryPattern:
                return relativedelta(minutes=1)
            elif any(hour in directoryPattern for hour in ["%H"]):
                return relativedelta(hours=1)
            elif any(day in directoryPattern for day in ["%d", "%j"]):
                return relativedelta(days=1)
            elif any(month in directoryPattern for month in ["%b", "%B", "%m"]):
                return relativedelta(months=1)
            elif any(year in directoryPattern for year in ["%Y", "%y"]):
                return relativedelta(years=1)
            else:
                return None
        except Exception:
            raise


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
    # This function deliberately does NOT use astropy.time because it is not
    # needed, and the performance overheads in dealing with astropy.time.Time
    # objects are large
    datetypes = ['year', 'month', 'day']
    timetypes = ['hour', 'minute', 'second', 'millisecond']
    dtlist = [int(exdict.get(d, 1)) for d in datetypes]
    dtlist.extend([int(exdict.get(t, 0)) for t in timetypes])
    startTime = datetime(*dtlist)
    tdelta = TIME_QUANTITIES['millisecond']
    if "second" in exdict:
        tdelta = TIME_QUANTITIES['second']
    elif "minute" in exdict:
        tdelta = TIME_QUANTITIES['minute']
    elif "hour" in exdict:
        tdelta = TIME_QUANTITIES['hour']
    elif "day" in exdict:
        tdelta = TIME_QUANTITIES['day']
    elif "month" in exdict:
        days_in_month = calendar.monthrange(int(exdict['year']), int(exdict['month']))[1]
        tdelta = days_in_month*TIME_QUANTITIES['day']
    elif "year" in exdict:
        if calendar.isleap(int(exdict['year'])):
            tdelta = 366*TIME_QUANTITIES['day']
        else:
            tdelta = 365*TIME_QUANTITIES['day']
    endTime = startTime + tdelta - TIME_QUANTITIES['millisecond']
    return TimeRange(startTime, endTime)
