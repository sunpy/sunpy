"""
This module provides a web scraper.
"""
import os
import re
from time import sleep
from ftplib import FTP
from datetime import datetime
from urllib.error import URLError, HTTPError
from urllib.parse import urlsplit
from urllib.request import urlopen

from bs4 import BeautifulSoup

from astropy.time import Time

from sunpy import log
from sunpy.extern.parse import parse
from sunpy.net.scraper_utils import date_floor, extract_timestep, get_timerange_from_exdict
from sunpy.time.timerange import TimeRange
from sunpy.util.decorators import deprecated_renamed_argument
from sunpy.util.exceptions import warn_user

__all__ = ['Scraper']


# Regular expressions to convert datetime format
TIME_CONVERSIONS = {
    '%Y': r'\d{4}', '%y': r'\d{2}',
    '%b': '[A-Z][a-z]{2}', '%B': r'\W', '%m': r'\d{2}',
    '%d': r'\d{2}', '%j': r'\d{3}',
    '%H': r'\d{2}', '%I': r'\d{2}',
    '%M': r'\d{2}',
    '%S': r'\d{2}', '%e': r'\d{3}', '%f': r'\d{6}'
}
# "parse" expressions to convert into datetime format
PARSE_TIME_CONVERSIONS = {
    "{year:4d}": "%Y", "{year:2d}": "%y",
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
DEPRECATED_MESSAGE = (
    "pattern has been replaced with the format keyword. "
    "This comes with a new syntax and there is a migration guide available at "
    "https://docs.sunpy.org/en/latest/topic_guide/scraper_migration.html."
)

class Scraper:
    """
    A scraper to scrap web data archives based on dates.

    Parameters
    ----------
    pattern : `str`
        A string containing the url with the date encoded as datetime formats,
        and any other parameter as ``kwargs`` as a string format.
        This can also be a uri to a local file patterns. Deprecated in favor of `format`. Default is `None`.
    regex : `bool`
        Set to `True` if parts of the pattern uses regexp symbols.
        This only works for the filename part of the pattern rather than the full url.
        Be careful that periods ``.`` matches any character and therefore it's better to escape them.
        If regexp is used, other ``kwargs`` are ignored and string replacement is
        not possible. Default is `False`.
        Deprecated in favor of `format`. Default is `False`.
    format : `str`
        A string containing the url with the date and other information to be
        extracted encoded as ``parse`` formats, and any other ``kwargs`` parameters
        as a string format, the former represented using double curly-brackets
        to differentiate from the latter.
        The accepted parse representations for datetime values are as given in ``PARSE_TIME_CONVERSIONS``.
        This can also be a uri to a local file patterns. Default is `None`.
    kwargs : `dict`
        A dictionary containing the values to be replaced in the pattern.
        Will be ignored if ``regex`` is `True`.

    Attributes
    ----------
    pattern : `str`
        The pattern with the parse format.
    datetime_pattern : `str`
        The parse pattern in the datetime format.
    now : `datetime.datetime`
        The pattern with the actual date.
        This is not checking if there is an existent file, but just how the ``pattern`` looks with the current time.

    Examples
    --------
    >>> from sunpy.net import Scraper
    >>>
    >>> pattern = ('http://proba2.oma.be/{instrument}/data/bsd/{{year:4d}}/{{month:2d}}/{{day:2d}}/'
    ...            '{instrument}_lv1_{{year:4d}}{{month:2d}}{{day:2d}}_{{hour:2d}}{{month:2d}}{{second:2d}}.fits')
    >>> swap = Scraper(format=pattern, instrument='swap')
    >>>
    >>> print(swap.pattern)
    http://proba2.oma.be/swap/data/bsd/{year:4d}/{month:2d}/{day:2d}/swap_lv1_{year:4d}{month:2d}{day:2d}_{hour:2d}{month:2d}{second:2d}.fits
    >>>
    >>> print(swap.datetime_pattern)
    http://proba2.oma.be/swap/data/bsd/%Y/%m/%d/swap_lv1_%Y%m%d_%H%m%S.fits
    >>>
    >>> print(swap.now)  # doctest: +SKIP
    http://proba2.oma.be/swap/data/bsd/2022/12/21/swap_lv1_20221221_112433.fits
    """
    @deprecated_renamed_argument("pattern", None, since="6.1", message=DEPRECATED_MESSAGE)
    @deprecated_renamed_argument("regex", None, since="6.1", message=DEPRECATED_MESSAGE)
    def __init__(self, pattern=None, regex=False, *, format=None, **kwargs):
        if pattern is not None and format is None:
            self.use_old_format = True
            self.pattern = pattern
        elif pattern is None and format is not None:
            self.use_old_format = False
            self.pattern = format
        else:
            raise ValueError("Either 'pattern' or 'format' must be provided.")
        if self.use_old_format:
            if regex and kwargs:
                warn_user('regexp being used, the extra arguments passed are being ignored')
            self.pattern = pattern.format(**kwargs) if kwargs and not regex else self.pattern
            self.domain = f"{urlsplit(self.pattern).scheme}://{urlsplit(self.pattern).netloc}/"
            milliseconds = re.search(r'\%e', self.pattern)
            if not milliseconds:
                self.now = datetime.now().strftime(self.pattern)
            else:
                now = datetime.now()
                milliseconds_ = int(now.microsecond / 1000.)
                self.now = now.strftime(f'{self.pattern[0:milliseconds.start()]}{milliseconds_:03d}{self.pattern[milliseconds.end():]}')
        else:
            pattern = format.format(**kwargs)
            datetime_pattern = pattern
            for k, v in PARSE_TIME_CONVERSIONS.items():
                if k in datetime_pattern:
                    datetime_pattern = datetime_pattern.replace(k,v)
            self.datetime_pattern = datetime_pattern
            # so that they don't conflict later on, either one can help in extracting year
            if "year:4d" in pattern and "year:2d" in pattern:
                pattern = pattern.replace("year:2d", ":2d")
            self.pattern = pattern
            self.domain = f"{urlsplit(self.pattern).scheme}://{urlsplit(self.pattern).netloc}/"
            milliseconds = re.search(r'\%e', self.datetime_pattern)
            if not milliseconds:
                self.now = datetime.now().strftime(self.datetime_pattern)
            else:
                now = datetime.now()
                milliseconds_ = int(now.microsecond / 1000.)
                self.now = now.strftime(f'{self.datetime_pattern[0:milliseconds.start()]}{milliseconds_:03d}{self.datetime_pattern[milliseconds.end():]}')

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
        if self.use_old_format:
            return date.strftime(self.pattern) == filepath
        else:
            return parse(date.strftime(self.datetime_pattern), filepath) is not None

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
        if self.use_old_format:
            if '/' in self.pattern:
                directorypattern = '/'.join(self.pattern.split('/')[:-1]) + '/'
        else:
            if '/' in self.datetime_pattern:
                directorypattern = '/'.join(self.datetime_pattern.split('/')[:-1]) + '/'
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
        >>> swap = Scraper(format=pattern, instrument='swap')
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
        >>> lyra = Scraper(format=pattern)
        >>> timerange = TimeRange('2023-03-06T00:00:00','2023-03-06T00:10:00')
        >>> print(swap.filelist(timerange)) # doctest: +REMOTE_DATA
        ['http://proba2.oma.be/swap/data/bsd/2023/03/06/swap_lv1_20230306_000128.fits',
        'http://proba2.oma.be/swap/data/bsd/2023/03/06/swap_lv1_20230306_000318.fits',
        'http://proba2.oma.be/swap/data/bsd/2023/03/06/swap_lv1_20230306_000508.fits',
        'http://proba2.oma.be/swap/data/bsd/2023/03/06/swap_lv1_20230306_000658.fits',
        'http://proba2.oma.be/swap/data/bsd/2023/03/06/swap_lv1_20230306_000848.fits']

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
        elif urlsplit(directories[0]).scheme in ["http", "https"]:
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
                    if self._url_follows_pattern(fullpath):
                        if self._check_timerange(fullpath, timerange):
                            filesurls.append(fullpath)

        filesurls = ['ftp://' + f"{urlsplit(url).netloc}{urlsplit(url).path}"
                     for url in filesurls]

        return filesurls

    def _localfilelist(self, timerange):
        """
        Goes over locally stored archives to return list of files in the given timerange.
        """
        if self.use_old_format:
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
                    if self._url_follows_pattern(fullpath):
                        if self._check_timerange(fullpath, timerange):
                            filepaths.append(fullpath)
            filepaths = [prefix + path for path in filepaths]
            self.pattern = pattern
            return filepaths
        else:
            pattern, datetime_pattern = self.pattern, self.datetime_pattern
            pattern_temp, datetime_pattern_temp = pattern.replace('file://', ''), datetime_pattern.replace('file://', '')
            if os.name == 'nt':
                pattern_temp = pattern_temp.replace('\\', '/')
                datetime_pattern_temp = datetime_pattern_temp.replace('\\', '/')
                prefix = 'file:///'
            else:
                prefix = 'file://'
            # Change pattern variables class-wide
            self.pattern, self.datetime_pattern = pattern_temp, datetime_pattern_temp
            directories = self.range(timerange)
            filepaths = list()
            for directory in directories:
                for file_i in os.listdir(directory):
                    fullpath = directory + file_i
                    if self._url_follows_pattern(fullpath):
                        if self._check_timerange(fullpath, timerange):
                            filepaths.append(fullpath)
            filepaths = [prefix + path for path in filepaths]
            # Set them back to their original values
            self.pattern, self.datetime_pattern = pattern, datetime_pattern
            return filepaths

    def _httpfilelist(self, timerange):
        """
        Goes over http archives hosted on the web, to return list of files in the given timerange.
        """
        directories = self.range(timerange)
        filesurls = list()
        retry_counts = {}
        while directories:
            directory = directories.pop(0)
            try:
                opn = urlopen(directory)
                try:
                    soup = BeautifulSoup(opn, "html.parser")
                    for link in soup.find_all("a"):
                        href = link.get("href")
                        if href is not None:
                            if href[0] == '/':
                                fullpath = self.domain + href[1:]
                            else:
                                fullpath = directory + href
                            if self._url_follows_pattern(fullpath):
                                if self._check_timerange(fullpath, timerange):
                                    filesurls.append(fullpath)
                finally:
                    opn.close()
            except HTTPError as http_err:
                # Ignore missing directories (issue #2684).
                if http_err.code == 404:
                    log.debug(f"Directory {directory} not found.")
                    continue
                if http_err.code in [400, 403]:
                    log.debug(f"Got error {http_err.code} while scraping {directory} : {http_err.reason}")
                    raise
                if http_err.code in [429, 504]:
                    # See if the server has told us how long to back off for
                    retry_after = http_err.hdrs.get('Retry-After', 2)
                    try:
                        # Ensure that we can parse the header as an int in sec
                        retry_after = int(retry_after)
                    except Exception as e:
                        log.debug(f"Converting retry_after failed: {e}")
                        retry_after = 2
                    log.debug(
                        f"Got {http_err.code} while scraping {directory}, waiting for {retry_after} seconds before retrying."
                    )
                    sleep(retry_after)
                    if retry_counts.get(directory, 0) > 4:
                     log.debug(f"Exceeded maximum retry limit for {directory}")
                     raise
                    retry_counts[directory] = retry_counts.get(directory, 0) + 1
                    directories.insert(0, directory)
                    continue
            except URLError as url_err:
               log.debug(f"Failed to parse content from {directory}: {url_err}")
               raise
            except Exception:
                raise
        return filesurls

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
        if self.use_old_format:
            if hasattr(self, 'extractor'):
                exdict = parse(self.extractor, url).named
                tr = get_timerange_from_exdict(exdict)
                return tr.intersects(timerange)
            else:
                datehref = self._extract_date(url).to_datetime()
                timestep = extract_timestep(self.pattern)
                tr = TimeRange(datehref, datehref + timestep)
                return tr.intersects(timerange)
        else:
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

    def _url_follows_pattern(self, url):
        """
        Check whether the url provided follows the pattern.
        """
        if self.use_old_format:
            pattern = self.pattern
            for k, v in TIME_CONVERSIONS.items():
                pattern = pattern.replace(k, v)
            matches = re.match(pattern, url)
            if matches:
                return matches.end() == matches.endpos
            return False
        else:
            return parse(self.pattern, url)


    def _extract_date(self, url):
        """
        Extracts the date from a particular url following the pattern.
        """
        # Remove the user and password if present
        url = url.replace("anonymous:data@sunpy.org@", "")

        def url_to_list(txt):
            # Substitutes '.' and '_' for '/'.
            return re.sub(r'\.|_', '/', txt).split('/')

        # Create a list of all the blocks in times
        # assuming they are all separated with either '.', '_' or '/'
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

    def _extract_files_meta(self, timerange, extractor=None, matcher=None):
        """
        Returns metadata information contained in URLs.

        Parameters
        ----------
        timerange : `~sunpy.time.TimeRange`
            Time interval where to find the directories for a given pattern.
        extractor : `str`
            Extractor to extract metadata from URLs if the old format for pattern is used.
        matcher : `dict`
            Dictionary to check if extracted metadata is valid.

        Returns
        -------
        `list` of `dict`
            List of metadata info for all URLs.
        """
        if self.use_old_format:
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
        else:
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
                            if match := matcher.get(k):
                                if str(metadict[k]) not in match:
                                    append = False
                                    break
                    if append:
                        metalist.append(metadict)
            return metalist
