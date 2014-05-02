import datetime

class Scrapper:
    def __init__(self, pattern, **kwargs):
        self.pattern = pattern.format(**kwargs)
        pass
    def matches(self, filepath, date):
        return date.strftime(self.pattern) == filepath
    def range(self, starttime, endtime):
        #find directory structure - without file names
        directorypattern = self.pattern[:-self.pattern[::-1].find('/')]
        #TODO what if there's not slashes?
        rangedelta = (endtime - starttime)
        timestep = self._smallerPattern(directorypattern)
        # Number of elements in the time range (including end)
        TotalTimeElements = int(rangedelta.total_seconds()/timestep.total_seconds()) + 1
        directories = [(starttime + n * timestep).strftime(directorypattern)
                       for n in range(TotalTimeElements)]
        return directories
    def _smallerPattern(self, directoryPattern):
        try: 
            if "%S" in directoryPattern:
                return datetime.timedelta(seconds=1)
            if "%M" in directoryPattern:
                return datetime.timedelta(minutes=1)
            if any(hour in directoryPattern for hour in ["%H", "%I"]):
                return datetime.timedelta(hours=1)
            if any(day in directoryPattern for day in ["%d", "%j"]):
                return datetime.timedelta(days=1)
            if any(month in directoryPattern for month in ["%b","%B","%m"]):
                return datetime.timedelta(days=31)
            if any(year in directoryPattern for year in ["%Y", "%y"]):
                return datetime.timedelta(days=365)
        except:
            raise 
