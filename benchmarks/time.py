from sunpy.time import is_time, parse_time


def time_is_time():
    """
    Benchmark the time taken to check if a string is a valid time.

    Example:
        is_time('1995-12-31 23:59:60')
    """
    is_time('1995-12-31 23:59:60')


def time_parse_time():
    """
    Benchmark the time taken to parse a time string into a SunPy-compatible time object.

    Example:
        parse_time('1995-12-31 23:59:60')
    """
    parse_time('1995-12-31 23:59:60')


def mem_parse_time():
    """
    Benchmark the memory usage when parsing a time string.

    Example:
        t = parse_time('1995-12-31 23:59:60')
        return t
    """
    t = parse_time('1995-12-31 23:59:60')
    return t


def peakmem_parse_time():
    """
    Benchmark the peak memory usage when parsing a time string.

    Example:
        parse_time('1995-12-31 23:59:60')
    """
    parse_time('1995-12-31 23:59:60')
