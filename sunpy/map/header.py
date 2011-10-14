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
        if isinstance(args[0], basestring):
            # filepath
            from sunpy.map.io import read_header
            #dict.__init__(self, dict(read_header(args[0])), **kwargs)
            dict.__init__(self, read_header(args[0]), **kwargs)
        else:
            # dictionary
            dict.__init__(self, *args, **kwargs)
            
    def copy(self):
        """Overide copy operator"""
        return type(self)(dict.copy(self))
        
    def __getitem__(self, key):
        """Overide [] indexing"""
        return dict.__getitem__(self, key.upper())
    
    def get(self, key, default=None):
        """Overide .get() indexing"""
        return dict.get(self, key.upper(), default)
    
if __name__ == "__main__":
    import sunpy
    aia = sunpy.Map(sunpy.AIA_171_IMAGE)
    print(aia[0:500,0:500])
        