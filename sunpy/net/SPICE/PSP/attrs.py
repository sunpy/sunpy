from sunpy.net.attr import SimpleAttr
from sunpy.time import parse_time


class Kernel_type(SimpleAttr):
    """
    kernel type
    """

class Analysis(SimpleAttr):
    """
    to get analysis frame kernel
    """
    def __init__(self,value):
        if not isinstance(value,bool):
            raise ValueError("given value must be boolean")
        
class Time(SimpleAttr):
    """
    Time attribute for kernel
    """

    def __init__(self,start,end = None):
        self.start = parse_time(start)
        self.end = parse_time(end) if end is not None else None

class Version(SimpleAttr):
    """
    version number for kernels
    """

class Numupdates(SimpleAttr):
    """
    Number of times file has been updated since launch
    """