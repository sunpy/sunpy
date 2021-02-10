
from sunpy.time import is_time


class TimeSuite:
    """
    An example benchmark that times the performance of sunpy.is_time()
    function
    """
    def time_is_time(self):
        t = is_time('1995-12-31 23:59:60')


class MemSuite:
    def mem_is_time(self):
        t = is_time('1995-11-31 23:59:60')
        return t
