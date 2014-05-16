import datetime
import urllib2
from bs4 import BeautifulSoup
import re

class Scrapper:
    def __init__(self, pattern, **kwargs):
        self.pattern = pattern.format(**kwargs)
        self.now = datetime.datetime.now().strftime(self.pattern)

    def matches(self, filepath, date):
        return date.strftime(self.pattern) == filepath

    def range(self, starttime, endtime):
        '''
        Gets the directories for a certain range of time.
        '''
        #find directory structure - without file names
        directorypattern = self.pattern[:-self.pattern[::-1].find('/')]
        #TODO what if there's not slashes?
        rangedelta = (endtime - starttime)
        timestep = self._smallerPattern(directorypattern)
        if timestep is None:
            return [directorypattern]
        else:
            # Number of elements in the time range (including end)
            TotalTimeElements = int(round(rangedelta.total_seconds()/timestep.total_seconds())) + 1
            directories = [(starttime + n * timestep).strftime(directorypattern)
                        for n in range(TotalTimeElements)] #todo if date <= endate
            return directories

    def URL_followsPattern(self, url):
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

    def extractDateURL(self, url):
        pattern_list = self.pattern.replace('.', '/').replace('_', '/').split('/')
        url_list = url.replace('.', '/').replace('_', '/').split('/')

        time_order = ['%Y', '%y', '%b', '%B', '%m', '%d', '%j', '%H', '%I', '%M', '%S']
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
        return datetime.datetime.strptime(' '.join(final_date), ' '.join(final_pattern))

    def filelist(self, starttime, endtime):
        directories = self.range(starttime, endtime)
        filesurls = []
        for directory in directories:
            opn = urllib2.urlopen(directory)
            try:
                soup = BeautifulSoup(opn)
                for link in soup.find_all("a"):
                    href = link.get("href")
                    if href.endswith(self.pattern.split('.')[-1]):
                        fullpath = directory + href
                        if self.URL_followsPattern(fullpath):
                            datehref = self.extractDateURL(fullpath)
                            if (datehref >= starttime and
                                datehref <= endtime):
                                filesurls.append(fullpath)
            finally:
                opn.close()
        return filesurls

    def _smallerPattern(self, directoryPattern):
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

