# Author: Simon Liedtke <liedtke.simon@googlemail.com>
#
# This module was developed with funding provided by
# the Google Summer of Code (2013).

from __future__ import absolute_import

from abc import ABCMeta, abstractmethod, abstractproperty
from collections import MutableMapping, OrderedDict, Counter

from sunpy.extern import six

__all__ = ['BaseCache', 'LRUCache', 'LFUCache']


@six.add_metaclass(ABCMeta)
class BaseCache(object):
    """
    BaseCache is a class that saves and operates on an OrderedDict. It has a
    certain capacity, stored in the attribute `maxsize`. Whether this
    capacity is reached, can be checked by using the boolean property
    `is_full`. To implement a custom cache, inherit from this class and
    override the methods ``__getitem__`` and ``__setitem__``.
    Call the method `sunpy.database.caching.BaseCache.callback` as soon
    as an item from the cache is removed.
    """

    def __init__(self, maxsize=float('inf')):
        self.maxsize = maxsize
        self._dict = OrderedDict()

    def get(self, key, default=None):  # pragma: no cover
        """Return the corresponding value to `key` if `key` is in the cache,
        `default` otherwise. This method has no side-effects, multiple calls
        with the same cache and the same passed key must always return the same
        value.

        """
        try:
            return self._dict[key]
        except KeyError:
            return default

    @abstractmethod
    def __getitem__(self, key):
        """abstract method: this method must be overwritten by inheriting
        subclasses. It defines what happens if an item from the cache is
        attempted to be accessed.

        """
        return  # pragma: no cover

    @abstractmethod
    def __setitem__(self, key, value):
        """abstract method: this method must be overwritten by inheriting
        subclasses. It defines what happens if a new value should be assigned
        to the given key. If the given key does already exist in the cache or
        not must be checked by the person who implements this method.
        """

    @abstractproperty
    def to_be_removed(self):
        """The item that will be removed on the next
        :meth:`sunpy.database.caching.BaseCache.remove` call.

        """

    @abstractmethod
    def remove(self):
        """Call this method to manually remove one item from the cache. Which
        item is removed, depends on the implementation of the cache. After the
        item has been removed, the callback method is called.

        """

    def callback(self, key, value):
        """This method should be called (by convention) if an item is removed
        from the cache because it is full. The passed key and value are the
        ones that are removed. By default this method does nothing, but it
        can be customized in a custom cache that inherits from this base class.

        """

    @property
    def is_full(self):
        """True if the number of items in the cache equals :attr:`maxsize`,
        False otherwise.

        """
        return len(self._dict) == self.maxsize

    def __delitem__(self, key):
        self._dict.__delitem__(key)

    def __contains__(self, key):
        return key in self._dict.keys()

    def __len__(self):
        return len(self._dict)

    def __iter__(self):
        for key in self._dict.__iter__():
            yield key

    def __reversed__(self):  # pragma: no cover
        for key in self._dict.__reversed__():
            yield key

    def clear(self):  # pragma: no cover
        return self._dict.clear()

    def keys(self):  # pragma: no cover
        return list(self._dict.keys())

    def values(self):  # pragma: no cover
        return list(self._dict.values())

    def items(self):  # pragma: no cover
        return list(self._dict.items())

    def iterkeys(self):  # pragma: no cover
        return iter(self._dict.keys())

    def itervalues(self):  # pragma: no cover
        for value in self._dict.values():
            yield value

    def iteritems(self):  # pragma: no cover
        for key, value in six.iteritems(self._dict):
            yield key, value

    def update(self, *args, **kwds):  # pragma: no cover
        self._dict.update(*args, **kwds)

    def pop(self, key, default=MutableMapping._MutableMapping__marker):  # pragma: no cover
        return self._dict.pop(key, default)

    def setdefault(self, key, default=None):  # pragma: no cover
        return self._dict.setdefault(key, default)

    def popitem(self, last=True):  # pragma: no cover
        return self._dict.popitem(last)

    def __reduce__(self):  # pragma: no cover
        return self._dict.__reduce__()

    def copy(self):  # pragma: no cover
        return self._dict.copy()

    def __eq__(self, other):  # pragma: no cover
        return self._dict.__eq__(other)

    def __ne__(self, other):  # pragma: no cover
        return self._dict.__ne__(other)

    def viewkeys(self):  # pragma: no cover
        return self._dict.keys()

    def viewvalues(self):  # pragma: no cover
        return self._dict.values()

    def viewitems(self):  # pragma: no cover
        return self._dict.items()

    @classmethod
    def fromkeys(cls, iterable, value=None):  # pragma: no cover
        return OrderedDict.fromkeys(iterable, value)

    def __repr__(self):  # pragma: no cover
        return '{0}({1!r})'.format(self.__class__.__name__, dict(self._dict))


class LRUCache(BaseCache):
    """
    LRUCache
    """
    @property
    def to_be_removed(self):
        """Return the least recently used key and its corresponding value as a
        tuple.

        """
        return six.next(self.iteritems())

    def remove(self):
        """Remove the least recently used item."""
        self.callback(*self.popitem(last=False))

    def __getitem__(self, key):
        """Returns the value which is associated to the given key and put it
        with its associated value to the end of this cache.

        Raises
        ------
        KeyError
            If the key cannot be found in the cache.

        """
        if key in self:
            value = self._dict.__getitem__(key)
            del self[key]
            self._dict.__setitem__(key, value)
            return value
        raise KeyError

    def __setitem__(self, key, value):
        """If the key does already exist in the cache, move it to the end of
        this cache. Otherwise, set a new value and put it to the end of this
        cache. If the cache is full, remove the least recently used item before
        inserting the new key-value pair.

        """
        if key in self:
            del self[key]
        if self.is_full:
            self.remove()
        self._dict.__setitem__(key, value)


class LFUCache(BaseCache):
    """
    LFUCache
    """
    def __init__(self, maxsize=float('inf')):
        self.usage_counter = Counter()
        BaseCache.__init__(self, maxsize)

    @property
    def to_be_removed(self):
        """Returns the key with the lowest times of access and its
        corresponding value as a tuple.

        """
        min_ = float('inf')
        lfu_key = None
        for k, v in six.iteritems(self.usage_counter):
            if v < min_:
                min_ = v
                lfu_key = k
        return lfu_key, self.get(lfu_key)

    def remove(self):
        """Remove the least frequently used item."""
        lfu_key, val = self.to_be_removed
        del self[lfu_key]
        del self.usage_counter[lfu_key]
        self.callback(lfu_key, val)

    def __getitem__(self, key):
        """Returns the value which is associated to the given key and
        increments the frequency counter of this key.

        Raises
        ------
        KeyError
            If the key cannot be found in the cache.

        """
        value = self._dict.__getitem__(key)
        self.usage_counter[key] += 1
        return value

    def __setitem__(self, key, value):
        """Increment the frequency counter of the given key if it is already
        present in the cache, otherwise set it to 1. If the cache is full,
        remove the least frequently used item before inserting the new
        key-value pair.

        """
        self.usage_counter[key] += 1
        if self.is_full:
            self.remove()
        self._dict.__setitem__(key, value)
