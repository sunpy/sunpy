from sunpy.net.attr import SimpleAttr
from sunpy.time import parse_time


#PSP SPECIFIC
class Analysis_fk(SimpleAttr):
    """
    to get analysis frame kernel
    """
    def __init__(self,value):
        self.value = value
        if not isinstance(value,bool):
            raise ValueError("given value must be boolean")

class Numupdates(SimpleAttr):
    """
    Number of times file has been updated since launch
    """

#Solo SPECIFIC
class Voem(SimpleAttr):
    """
    Voem : reference to the source OEM file version
    """

class Readme(SimpleAttr):
    def __init__(self,value):

        if not isinstance(value,bool):
            raise ValueError(f"value must be boolean not {type(value)}")
        self.value = value

class Sensor(SimpleAttr):
    """
    sensor for kernels
    """

#general
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
