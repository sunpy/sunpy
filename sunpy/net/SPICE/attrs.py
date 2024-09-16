from sunpy.net.attr import SimpleAttr
from sunpy.time import parse_time


class Kernel_type(SimpleAttr):
    """
    kernel type
    """
    def __init__(self, value):
        self.value = value
        if value is None:
            raise ValueError ("kernel type is required")

class Instrument(SimpleAttr):
    """
    Instrument for kernels
    """

class Link(SimpleAttr):
    """
    name of link for spice kernels
    """

class Version(SimpleAttr):
    """
    version number for kernels
    """

class Time(SimpleAttr):
    """
    Time attribute for kernel
    """
    def __init__(self,start,end = None):
        self.start = parse_time(start)
        self.end = parse_time(end) if end is not None else None

class Index(SimpleAttr):
    """
    index of link to be downloaded
    """
    def __init__(self,*value):
        self.value = value

class Mission(SimpleAttr):
    def __init__(self, *value):
        self.value = value
