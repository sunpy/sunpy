from sunpy.net.attr import SimpleAttr

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

class sensor(SimpleAttr):
    """
    sensor for kernels
    """
class Instrument(SimpleAttr):
    """
    Instrument for kernels
    """

class link(SimpleAttr):
    """
    kernel links
    """