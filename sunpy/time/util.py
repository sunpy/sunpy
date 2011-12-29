import datetime
from sunpy.util.util import anytim

class TimeRange:
    """A class to hold a time range"""
    def __init__(self, t1, t2):
        self.t1 = anytim(t1)
        self.t2 = anytim(t2)
        self.dt = self.t2 - self.t1
            
    def days(self):
        return self.dt.days
    
    def seconds(self):
        return self.dt.total_seconds()
    
    def minutes(self):
        return self.dt.total_seconds()/60.0