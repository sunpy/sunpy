from sunpy.net.attr import SimpleAttr


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
