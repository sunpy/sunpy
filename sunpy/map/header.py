"""
MapHeader is a generalized header class that deals with header parsing and
normalization.
"""
from __future__ import absolute_import

from collections import OrderedDict

__all__ = ['MapMeta']

class MapMeta(OrderedDict):
    """
    A class to hold meta data associated with a Map derivative.
    
    This class handles everything a lower case. This allows case insensitive 
    indexing.
    """
    def __init__(self, adict):
        """Creates a new MapHeader instance"""
        # Store all keys as upper-case to allow for case-insensitive indexing
        adict = dict((k.lower(), v) for k, v in adict.items())
        OrderedDict.__init__(self, adict)
        
    def __contains__(self, key):
        """Overide __contains__"""
        return dict.__contains__(self, key.lower())

    def __getitem__(self, key):
        """Overide [] indexing"""
        return dict.__getitem__(self, key.lower())

    def __setitem__(self, key, value):
        """Overide [] indexing"""
        return dict.__setitem__(self, key.lower(), value)
    
    def as_pyfits_header(self):
        """Returns a PyFITS header instance of the header"""
        cards = [pyfits.core.Card(k, v) for k, v in self.items()]
        return pyfits.core.Header(cards)

    def copy(self):
        """Overide copy operator"""
        return type(self)(dict.copy(self))

    def get(self, key, default=None):
        """Overide .get() indexing"""
        return dict.get(self, key.lower(), default)

    def has_key(self, key):
        """Overide .has_key() to perform case-insensitively"""
        return key.lower() in self

    def pop(self, key, default=None):
        """Overide .pop() to perform case-insensitively"""
        return dict.pop(self, key.lower(), default)

    def update(self, d2):
        """Overide .update() to perform case-insensitively"""
        return dict.update(self, dict((k.lower(), v) for k, v in d2.items()))

    def setdefault(self, key, default=None):
        """Overide .setdefault() to perform case-insensitively"""
        return dict.setdefault(self, key.lower(), default)