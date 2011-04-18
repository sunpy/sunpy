"""Map class creation functions

Author: `Keith Hughitt <keith.hughitt@nasa.gov>`_

"""
__author__ = "Keith Hughitt"
__email__ = "keith.hughitt@nasa.gov"

import pyfits
from BaseMap import BaseMap
from sources import *

def Map(input_):
    """Map class factory
    
    Attempts to determine the type of data associated with input and returns
    an instance of either the generic BaseMap class or a subclass of BaseMap
    such as AIAMap, EUVIMap, etc.
    
    Parameters
    ----------
    input_ : filepath, data array
        The data source used to create the map object. This can be either a
        filepath to an image, a 2d list, or an ndarray.
        
    Returns
    -------
    out : BaseMap
        Returns a BaseMap or BaseMap subclass instance
    
    See Also
    --------
    http://stackoverflow.com/questions/456672/class-factory-in-python
    """
    if isinstance(input_, str):
        try:
            fits = pyfits.open(input_)
            data = fits[0].data
            header = fits[0].header
        except IOError:
            sys.exit("Unable to read the file %s" % input_)

        for cls in BaseMap.__subclasses__():
            if cls.is_datasource_for(header):
                return cls(data, header)
        raise UnrecognizedDataSouce

    else:
        return BaseMap(input_)
    
class UnrecognizedDataSouce(ValueError):
    """Exception to raise when an unknown datasource is encountered"""
    pass
