
class Scrapper:
    def __init__(self, pattern, **kwargs):
        self.pattern = pattern.format(**kwargs)
        pass
    def matches(self, filepath, date):
        return date.strftime(self.pattern) == filepath
    def range(self, startime, endtime):
        #find directory structure
        directorypattern = self.pattern[:-self.pattern[::-1].find('/')]
        #directories = [d.strftime(directorypattern) for d in ]
        # we may need to remove duplicates
        return True
