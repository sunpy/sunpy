"""
This module provides a generalized dictionary class that deals with header
parsing, normalization, and maintaining coherence between keys and keycomments.
"""
from collections import OrderedDict, namedtuple
from collections.abc import Mapping

__all__ = ['MetaDict']

ModifiedItem = namedtuple('ModifiedItem', ['original', 'current'])
ModifiedItem.__repr__ = lambda t: f"(original={t.original}, current={t.current})"


class MetaDict(OrderedDict):
    """
    A class to hold metadata associated with a `sunpy.map.Map` derivative.

    This class handles everything in lower case. This allows case
    insensitive indexing.

    If the key 'keycomments' exists, its value must be a dictionary mapping
    keys in the :class:`MetaDict` to their comments. The casing of keys in the
    keycomments dictionary is not significant. If a key is removed from the
    :class:`MetaDict`, it will also be removed from the keycomments dictionary.
    Additionally, any extraneous keycomments will be removed when the
    :class:`MetaDict` is instantiated.

    Parameters
    ----------
    save_original : bool, optional
        If `True` (the default), a copy of the original metadata will be saved to the
        `original_meta` property. Note that this keyword argument is
        required as an implementation detail to avoid recursion when saving
        the original contents.
    """

    def __init__(self, *args, save_original=True):
        # Store all keys as lower-case to allow for case-insensitive indexing
        # OrderedDict can be instantiated from a list of lists or a tuple of tuples
        tags = dict()
        if args:
            args = list(args)
            adict = args[0]

            if isinstance(adict, MetaDict):
                self._original_meta = adict._original_meta

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

        if save_original and not hasattr(self, '_original_meta'):
            self._original_meta = MetaDict(*args, save_original=False)

    def __str__(self):
        return '\n'.join([f"('{key}': '{item}')" for key, item in self.items()])

    def __repr__(self):
        return f"{self.__class__.__name__}([{self}])"

    # Deliberately a property to prevent external modification
    @property
    def original_meta(self):
        """
        A copy of the dictionary from when it was initialised.
        """
        return self._original_meta

    @property
    def added_items(self):
        """
        Items added since initialisation.
        """
        return {k: self[k] for k in set(self) - set(self.original_meta)}

    @property
    def removed_items(self):
        """
        Items removed since initialisation.
        """
        return {k: self.original_meta[k] for k in set(self.original_meta) - set(self)}

    @property
    def modified_items(self):
        """
        Items modified since initialisation.

        This is a mapping from ``key`` to ``[original_value, current_value]``.
        """
        keys = set(self).intersection(set(self.original_meta))
        return {k: ModifiedItem(self.original_meta[k], self[k]) for k in keys if
                self[k] != self.original_meta[k]}

    def copy(self):
        """
        Shallow copies this instance.
        """
        copied = super().copy()
        # By default the above line will overwrite original_meta, so manually re-instate it
        copied._original_meta = self.original_meta
        return copied

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
                f"the following type: {type(keycomments)!r}")

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

    def popitem(self, last):
        key, value = super().popitem(last)
        self._prune_keycomments()
        return key, value

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
