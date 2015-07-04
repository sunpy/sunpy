import datetime
import re

import urllib2
from bs4 import BeautifulSoup

__all__ = ['Scraper']

class Scraper(object):
    """
    A Scraper to scrap web data archives based on dates.

    Parameters
    ----------
    pattern : string
        A string containing the url with the date encoded as
        datetime formats, and any other parameter as kwargs
        as string format.

    Attributes
    ----------
    pattern : string
        A converted string with the kwargs.
    now : datetime.datetime
        The pattern with the actual date.

    Examples
    --------
    >>> from sunpy.util.scraper import Scraper
    >>> solmon_pattern = ('http://solarmonitor.org/data/'
                          '%Y/%m/%d/fits/{instrument}/'
                          '{instrument}_{wave:05d}_fd_%Y%m%d_%H%M%S.fts.gz')
    >>> solmon = Scraper(solmon_pattern, instrument = 'swap', wave = 174)
    >>> print(solmon.pattern)
    http://solarmonitor.org/data/%Y/%m/%d/fits/swap/swap_00174_fd_%Y%m%d_%H%M%S.fts.gz
    >>> print(solmon.now)
    http://solarmonitor.org/data/2012/01/25/fits/swap/swap_00174_fd_20120125_173301.fts.gz

    Notes
    -----
    The now attribute does not return an existent file, but just how the
    pattern looks with the actual time.
    """
    def __init__(self, pattern, **kwargs):
        self.pattern = pattern.format(**kwargs)
        self.now = datetime.datetime.now().strftime(self.pattern)

    def matches(self, filepath, date):
        return date.strftime(self.pattern) == filepath

    def range(self, timerange):
        """
        Gets the directories for a certain range of time
        (i.e. using `~sunpy.time.TimeRange`).

        Parameters
        ----------

        timerange : `~sunpy.time.TimeRange`
            Time interval where to find the directories for a given
            pattern.

        Returns
        -------

        directories : list of strings
            List of all the possible directories valid for the time
            range given. Notice that these directories may not exist
            in the archive.
        """
        #find directory structure - without file names
        directorypattern = self.pattern[:-self.pattern[::-1].find('/')]
        #TODO what if there's not slashes?
        rangedelta = timerange.dt # TODO check it's a timerange
        timestep = self._smallerPattern(directorypattern)
        if timestep is None:
            return [directorypattern]
        else:
            # Number of elements in the time range (including end)
            n_steps = rangedelta.total_seconds()/timestep.total_seconds()
            TotalTimeElements = int(round(n_steps)) + 1
            directories = [(timerange.start + n * timestep).strftime(directorypattern)
                        for n in range(TotalTimeElements)] #todo if date <= endate
            return directories

    def _URL_followsPattern(self, url):
        """Check whether the url provided follows the pattern"""
        time_conversions = {'%Y': '\d{4}', '%y': '\d{2}',
                            '%b': '[A-Z]..', '%B': '\W', '%m': '\d{2}',
                            '%d': '\d{2}', '%j': '\d{3}',
                            '%H': '\d{2}', '%I': '\d{2}',
                            '%M': '\d{2}',
                            '%S': '\d{2}'}
        pattern = self.pattern
        for k,v in time_conversions.iteritems():
            pattern = pattern.replace(k, v)
        matches = re.match(pattern, url)
        if matches:
            return matches.end() == matches.endpos == len(self.now)
        return False

    def _extractDateURL(self, url):
        """Extracts the date from a particular url following the pattern"""
        url_to_list = lambda txt: re.sub(r'\.|_', '/', txt).split('/')
        pattern_list = url_to_list(self.pattern)
        url_list = url_to_list(url)

        time_order = ['%Y', '%y', '%b', '%B', '%m', '%d', '%j',
                      '%H', '%I', '%M', '%S']
        final_date = []
        final_pattern = []
        #Find in directory and filename
        for pattern_elem, url_elem in zip(pattern_list, url_list):
            time_formats = [x for x in time_order if x in pattern_elem]
            if len(time_formats) > 0:
                final_date.append(url_elem)
                final_pattern.append(pattern_elem)
                for time_bit in time_formats:
                    time_order.remove(time_bit)
        try:
            return datetime.datetime.strptime(' '.join(final_date), ' '.join(final_pattern))
        except:
            anomaly = [x for x in final_pattern if len(x) > 2]
            for x in anomaly:
                splits = [x[i:i+2] for i in range(0, len(x), 2)]
                index_remove = [final_pattern.index(i) for i in splits if i in final_pattern]
                for i in splits:
                    if i in final_pattern: final_pattern.remove(i)
                for i in index_remove: final_date.remove(final_date[i])
            return datetime.datetime.strptime(' '.join(final_date), ' '.join(final_pattern)) 

    def filelist(self, timerange):
        """
        Returns the list of existent files in the archive for the
        given time range.

        Parameters
        ----------

        timerange : `~sunpy.time.TimeRange`
            Time interval where to find the directories for a given
            pattern.

        Returns
        -------

        filesurls : list of strings
            List of all the files found between the time range given.

        Examples
        --------
        >>> from sunpy.time import TimeRange
        >>> timerange = TimeRange('2015-01-01','2015-01-01T16:00:00')
        >>> print(solmon.filelist(timerange))
        ['http://solarmonitor.org/data/2015/01/01/fits/swap/swap_00174_fd_20150101_025423.fts.gz']
        """
        directories = self.range(timerange)
        filesurls = []
        for directory in directories:
            try:
                opn = urllib2.urlopen(directory)
                try:
                    soup = BeautifulSoup(opn)
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
            except:
                pass
        return filesurls

    def _smallerPattern(self, directoryPattern):
        """Obtain the smaller time step for the given pattern"""
        try: 
            if "%S" in directoryPattern:
                return datetime.timedelta(seconds=1)
            elif "%M" in directoryPattern:
                return datetime.timedelta(minutes=1)
            elif any(hour in directoryPattern for hour in ["%H", "%I"]):
                return datetime.timedelta(hours=1)
            elif any(day in directoryPattern for day in ["%d", "%j"]):
                return datetime.timedelta(days=1)
            elif any(month in directoryPattern for month in ["%b","%B","%m"]):
                return datetime.timedelta(days=31)
            elif any(year in directoryPattern for year in ["%Y", "%y"]):
                return datetime.timedelta(days=365)
            else:
                return None
        except:
            raise

