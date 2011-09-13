from __future__ import absolute_import
"""
MapHeader is a generalized header class that deals with header parsing and
normalization.
"""

class MapHeader(dict):
    """
    MapHeader(header)
    
    A dictionary-like class for working with FITS, etc headers
    
    Parameters
    ----------
    header : pyfits.core.Header, dict
        Header tags associated with the data
        
    Attributes
    ----------
    
    """
    def __init__(self, *args, **kwargs):
        """Creates a new MapHeader instance"""
        dict.__init__(self, *args, **kwargs)
        
    def __getitem__(self, key):
        """Overide [] indexing"""
        return dict.__getitem__(self, key.upper())
    
    def get(self, key, default=None):
        """Overide .get() indexing"""
        return dict.get(self, key.upper(), default)
        