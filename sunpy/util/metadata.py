"""
MapHeader is a generalized header class that deals with header parsing and
normalization.
"""
from __future__ import absolute_import, division, print_function

from collections import OrderedDict

__all__ = ['MetaDict']


class MetaDict(OrderedDict):
    """
    A class to hold meta data associated with a Map derivative.

    This class handles everything in lower case. This allows case insensitive
    indexing.
    """
    def __init__(self, *args):
        """Creates a new MapHeader instance"""
        # Store all keys as upper-case to allow for case-insensitive indexing
        # OrderedDict can be instantiated from a list of lists or a tuple of tuples
        tags = dict()
        if args:
            args = list(args)
            adict = args[0]
            if isinstance(adict, list) or isinstance(adict, tuple):
                tags = OrderedDict((k.upper(), v) for k, v in adict)
            elif isinstance(adict, dict):
                tags = OrderedDict((k.upper(), v) for k, v in adict.items())
            else:
                raise TypeError("Can not create a MetaDict from this type input")
            args[0] = tags

        super(MetaDict, self).__init__(*args)

    def __contains__(self, key):
        """Override __contains__"""
        return OrderedDict.__contains__(self, key.lower())

    def __getitem__(self, key):
        """Override [] indexing"""
        return OrderedDict.__getitem__(self, key.lower())

    def __setitem__(self, key, value):
        """Override [] indexing"""
        return OrderedDict.__setitem__(self, key.lower(), value)

    def get(self, key, default=None):
        """Override .get() indexing"""
        return OrderedDict.get(self, key.lower(), default)

    def has_key(self, key):
        """Override .has_key() to perform case-insensitively"""
        return key.lower() in self

    def pop(self, key, default=None):
        """Override .pop() to perform case-insensitively"""
        return OrderedDict.pop(self, key.lower(), default)

    def update(self, d2):
        """Override .update() to perform case-insensitively"""
        return OrderedDict.update(self, OrderedDict((k.lower(), v) for k, v in d2.items()))

    def setdefault(self, key, default=None):
        """Override .setdefault() to perform case-insensitively"""
        return OrderedDict.setdefault(self, key.lower(), default)
