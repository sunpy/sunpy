from sunpy.net.attr import SimpleAttr


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
