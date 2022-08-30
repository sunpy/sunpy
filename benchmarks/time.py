from sunpy.time import is_time, parse_time


def time_is_time(self):
    is_time('1995-12-31 23:59:60')


def time_parse_time(self):
    parse_time('1995-12-31 23:59:60')


def mem_parse_time(self):
    t = parse_time('1995-12-31 23:59:60')
    return t


def peakmem_parse_time(self):
    t = parse_time('1995-12-31 23:59:60')
