"""Definition for a Python Map Object

Author: Keith Hughitt <keith.hughitt@nasa.gov>

See: http://docs.scipy.org/doc/numpy/reference/arrays.classes.html
"""

import numpy as np
import pyfits
from datetime import datetime

class Map(np.ndarray):
    """A spatially-aware data array"""
    def __new__(cls, input_):
        """Creates a new Map instance
        
        Parameters
        ----------
        input_ : {filepath, data array}
            The data source used to create the map object
        """
        if isinstance(input_, str):
            try:
                fits = pyfits.open(input_)
            except IOError:
                sys.exit("Unable to read the file %s" % input_)

            #np.ndarray.__new__(self, fits[0].data)
            obj = np.asarray(fits[0].data).view(cls)
            
            obj.header = fits[0].header
            
            norm_header = parse_header(obj.header)
            
            obj.date = norm_header['date']
            obj.det = norm_header['det']
            obj.inst = norm_header['inst']
            obj.meas = norm_header['meas']
            obj.obs = norm_header['obs']
            obj.name = norm_header['name']
            obj.r_sun = norm_header['r_sun']

        elif isinstance(input_, list):
            obj = np.asarray(input_).view(cls)
        elif isinstance(input_, np.ndarray):
            obj = input_
            
        return obj
            
    def __array_finalize__(self, obj):
        """Finishes instantiation of the new map object"""
        if obj is None: return

#
# External functions - These should be moved to a separate file or files
# once an appropriate location is determined
#         
def parse_header(header):
    """Parses a FITS, etc image header
    
    Attempts to detect the type of data (e.g. AIA) based on values in the 
    file header and returns a mapped dictionary of of some important values.

    Parameters
    ----------
    header : {dict} A dictionary container the header keywords from the file 
        being read in
    """
    # AIA
    if header['telescop'] == 'SDO/AIA':
        return get_norm_header_tags(header, "aia")

def get_norm_header_tags(header, type_):
    """Returns a normalized dictionary of header values
    
    A normalized mapping of important header values is created and returned.
    Not all of the header values are used, but instead only those that are
    required for the Map class to function are included. Note that some values
    may be cast to new types in the process.
    
    Parameters
    ----------
    header : {dict} A dictionary container the header keywords from the file 
        being read in
    type_: A short string describing the type of data being mapped
    
    Returns
    -------
    out: {dict} A new mapped dictionary of useful header values
    """
    date_fmt1 = "%Y-%m-%dT%H:%M:%S.%f"
    
    if type_ == "aia":
        return {
            "date": datetime.strptime(header['date-obs'], date_fmt1),
            "det": "AIA",
            "inst": "AIA",
            "meas": 171,
            "obs": "SDO",
            "name": "AIA %s" % header['wave_str'][0:3],
            "r_sun": header['r_sun']
        }
