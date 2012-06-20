"""
LightCurve is a generic LightCurve class from which all other LightCurve classes 
inherit from.
"""
from __future__ import absolute_import

#pylint: disable=E1101,E1121,W0404,W0613
__authors__ = ["Keith Hughitt"]
__email__ = "keith.hughitt@nasa.gov"

from sunpy.io import read_file

class LightCurve:
    """
    LightCurve(filepath)

    A generic light curve object.

    Parameters
    ----------
    filepath : string
        Location of file to read in

    Attributes
    ----------


    Examples
    --------

    See Also
    --------
    

    References
    ----------


    """
    def __init__(self, data, header):
        self.data = data
        self.header = header
        
    @classmethod
    def parse_file(cls, filepath):
        """Reads in a map file and returns a header and data array"""
        data, dict_header = read_file(filepath)

        return dict_header, data

    @classmethod
    def read(cls, filepath):
        """LightCurve class factory

        Attempts to determine the type of data associated with input and
        returns a LightCurve subclass instance. If the file format is not
        recognized a warning will be displayed.

        Parameters
        ----------
        filepath : string
            Path to input file (FITS, CSV, etc.)

        Returns
        -------
        out : LightCurve
            Returns a LightCurve instance.
        """
        header, data = cls.parse_file(filepath)

        if cls.__name__ is not "LightCurve":
            return cls(filepath)

        for cls in LightCurve.__subclasses__():
            if cls.is_datasource_for(header):
                return cls(data, header)
        
        # DISPLAY WARNING..
