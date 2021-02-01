
from sunpy.time import parse_time


class TimeSuite:
    """
    An example benchmark that times the performance of sunpy.parse_time()
    function
    """
    def time_parse_time(self):
        t = parse_time('1995-12-31 23:59:60')


class MemSuite:
    def mem_parse_time(self):
        t = parse_time('1995-12-31 23:59:60')
        return t


class PeakMemorySuite:
    def peakmem_parse_time(self):
        t = parse_time('1995-12-31 23:59:60')
