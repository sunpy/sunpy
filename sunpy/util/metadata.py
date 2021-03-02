"""
This module provides a generalized dictionary class that deals with header
parsing, normalization, and maintaining coherence between keys and keycomments.
"""
from collections import OrderedDict
from collections.abc import Mapping

__all__ = ['MetaDict']


class MetaDict(OrderedDict):
    """
    A class to hold metadata associated with a `sunpy.map.Map` derivative.

    This class handles everything in lower case. This allows case
    insensitive indexing.

    If the key 'keycomments' exists, its value must be a dictionary mapping
    keys in the `MetaDict` to their comments. The casing of keys in the
    keycomments dictionary is not significant. If a key is removed from the
    `MetaDict`, it will also be removed from the keycomments dictionary.
    Additionally, any extraneous keycomments will be removed when the
    `MetaDict` is instantiated.
    """

    def __init__(self, *args):
        """
        Creates a new MetaDict instance.
        """
        # Store all keys as lower-case to allow for case-insensitive indexing
        # OrderedDict can be instantiated from a list of lists or a tuple of tuples
        tags = dict()
        if args:
            args = list(args)
            adict = args[0]
            if isinstance(adict, list) or isinstance(adict, tuple):
                items = adict
            elif isinstance(adict, Mapping):
                # Cast to a dict here, since a Mapping doesn't have to define
                # a .items() method (but has enough methods to be converted to
                # a dict)
                items = dict(adict).items()
            else:
                raise TypeError(f"Can not create a MetaDict from this input of type {type(adict)}")

            self._check_str_keys(items)
            tags = OrderedDict((k.lower(), v) for k, v in items)
            args[0] = tags

        super().__init__(*args)

        # Use `copy=True` to avoid mutating the caller's keycomments
        # dictionary (if they provided one).
        self._prune_keycomments(copy=True)

    @staticmethod
    def _check_str_keys(items):
        bad_keys = []
        for k, v in items:
            if not isinstance(k, str):
                bad_keys.append(str(k))
        if len(bad_keys):
            raise ValueError('All MetaDict keys must be strings, '
                             'found the following non-compliant keys: ' +
                             ', '.join(bad_keys))

    def _prune_keycomments(self, copy=False):
        """
        Remove keycomments for keys that are not contained in the MetaDict.

        Parameters
        ----------
        copy : `bool`, optional
            Make a copy of the current keycomments dict before removing keys.
        """
        if 'keycomments' not in self:
            return

        keycomments = self['keycomments']

        if not isinstance(keycomments, dict):
            raise TypeError(
                "'keycomments' key must have a value of type `dict`. Found "
                "the following type: %r" % type(keycomments))

        if copy:
            keycomments = keycomments.copy()

        for key in list(keycomments.keys()):
            if key not in self:
                del keycomments[key]

        self['keycomments'] = keycomments

    def __contains__(self, key):
        """
        Override ``__contains__``.
        """
        return OrderedDict.__contains__(self, key.lower())

    def __getitem__(self, key):
        """
        Override ``[]`` indexing.
        """
        return OrderedDict.__getitem__(self, key.lower())

    def __setitem__(self, key, value):
        """
        Override ``[]`` indexing.
        """
        return OrderedDict.__setitem__(self, key.lower(), value)

    # Note: `OrderedDict.popitem()` does not need to be overridden to prune
    # keycomments because it calls `__delitem__` internally.
    def __delitem__(self, key):
        """
        Override ``del dict[key]`` key deletion.
        """
        OrderedDict.__delitem__(self, key.lower())
        self._prune_keycomments()

    def item_hash(self):
        """
        Create a hash based on the stored items.

        This relies on all the items themselves being hashable. For this reason
        the 'keycomments' item, which is a dict, is excluded from the hash.

        If creating the hash fails, returns `None`.
        """
        self_copy = self.copy()
        self_copy.pop('keycomments', None)
        try:
            return hash(frozenset(self_copy.items()))
        except Exception:
            return

    def get(self, key, default=None):
        """
        Override ``.get()`` indexing.
        """
        return OrderedDict.get(self, key.lower(), default)

    def has_key(self, key):
        """
        Override ``.has_key()`` to perform case-insensitively.
        """
        return key.lower() in self

    def pop(self, key, default=None):
        """
        Override ``.pop()`` to perform case-insensitively.
        """
        has_key = key in self
        result = OrderedDict.pop(self, key.lower(), default)
        if has_key:
            self._prune_keycomments()
        return result

    def update(self, d2):
        """
        Override ``.update()`` to perform case-insensitively.
        """
        return OrderedDict.update(self, OrderedDict((k.lower(), v) for k, v in d2.items()))

    def setdefault(self, key, default=None):
        """
        Override ``.setdefault()`` to perform case-insensitively.
        """
        return OrderedDict.setdefault(self, key.lower(), default)
