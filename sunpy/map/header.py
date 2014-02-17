"""
MapHeader is a generalized header class that deals with header parsing and
normalization.
"""
from __future__ import absolute_import

from sunpy.util.odict import OrderedDict

__all__ = ['MapMeta']

class MapMeta(OrderedDict):
    """
    A class to hold meta data associated with a Map derivative.
    
    This class handles everything a lower case. This allows case insensitive 
    indexing.
    """
    def __init__(self, adict, *args):
        """Creates a new MapHeader instance"""
        # Store all keys as upper-case to allow for case-insensitive indexing
        #OrderedDict can be instanciated from a list of lists or a tuple of tuples
        if isinstance(adict, list) or isinstance(adict, tuple):
            tags = dict((k.upper(), v) for k, v in adict)
        elif isinstance(adict, dict):
            tags = dict((k.upper(), v) for k, v in adict.items())
        else:
            raise TypeError("Can not create a MapMeta from this type input")
            
        super(MapMeta, self).__init__(tags, *args)

    def __contains__(self, key):
        """Overide __contains__"""
        return OrderedDict.__contains__(self, key.lower())

    def __getitem__(self, key):
        """Overide [] indexing"""
        return OrderedDict.__getitem__(self, key.lower())

    def __setitem__(self, key, value):
        """Overide [] indexing"""
        return OrderedDict.__setitem__(self, key.lower(), value)
        
    def get(self, key, default=None):
        """Overide .get() indexing"""
        return OrderedDict.get(self, key.lower(), default)

    def has_key(self, key):
        """Overide .has_key() to perform case-insensitively"""
        return key.lower() in self

    def pop(self, key, default=None):
        """Overide .pop() to perform case-insensitively"""
        return OrderedDict.pop(self, key.lower(), default)

    def update(self, d2):
        """Overide .update() to perform case-insensitively"""
        return OrderedDict.update(self, dict((k.lower(), v) for k, v in d2.items()))

    def setdefault(self, key, default=None):
        """Overide .setdefault() to perform case-insensitively"""
        return OrderedDict.setdefault(self, key.lower(), default)