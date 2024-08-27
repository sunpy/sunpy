from sunpy.net.attr import SimpleAttr
from sunpy.time import parse_time


class Kernel_type(SimpleAttr):
    """
    kernel type
    """
    def __init__(self, value):
        super().__init__(value)
        if value is None:
            raise ValueError ("kernel type is required")
        if value not in ["ck", "fk", "ik", "lsk", "pck", "sclk", "spk","mk"]:
            raise ValueError(f"Kernel type not recognized '{value}'")

class Sensor(SimpleAttr):
    """
    sensor for kernels
    """
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


class Readme(SimpleAttr):
    def __init__(self,value):

        if not isinstance(value,bool):
            raise ValueError(f"value must be boolean not {type(value)}")
        self.value = value

class Voem(SimpleAttr):
    """
    Voem : reference to the source OEM file version
    """
