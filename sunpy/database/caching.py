from __future__ import absolute_import

from abc import ABCMeta, abstractmethod, abstractproperty
from collections import OrderedDict, Counter

__all__ = ['BaseCache', 'LRUCache', 'LFUCache']


class BaseCache(object):
    """BaseCache is a class that saves and operates on an OrderedDict. It has a
    certain capacity, stored in the attribute :attr:`maxsize`. Whether this
    capacity is reached, can be checked by using the boolean property
    :attr:`is_full`. To implement a custom cache, inherit from this class and
    override the methods ``__getitem__`` and ``__setitem__``.
    Call the method :meth:`callback` as soon as an item from the cache is
    removed.

    Methods
    -------
    get
    remove
    callback

    """
    __metaclass__ = ABCMeta

    def __init__(self, maxsize=float('inf')):
        self.maxsize = maxsize
        self.dict = OrderedDict()

    def get(self, key, default=None):  # pragma: no cover
        """Return the corresponding value to `key` if `key` is in the cache,
        `default` otherwise. This method has no side-effects, multiple calls
        with the same cache and the same passed key must always return the same
        value.

        """
        try:
            return self.dict[key]
        except KeyError:
            return default

    @abstractmethod
    def __getitem__(self, key):
        return  # pragma: no cover

    @abstractmethod
    def __setitem__(self, key, value):
        return  # pragma: no cover

    @abstractproperty
    def to_be_removed(self):
        return  # pragma: no cover

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
        return len(self.dict) == self.maxsize

    def __delitem__(self, key):
        self.dict.__delitem__(key)

    def __contains__(self, key):
        return key in self.dict.keys()

    def __len__(self):
        return len(self.dict)

    def __iter__(self):
        for key in self.dict.__iter__():
            yield key

    def __reversed__(self):  # pragma: no cover
        for key in self.dict.__reversed__():
            yield key

    def clear(self):  # pragma: no cover
        return self.dict.clear()

    def keys(self):  # pragma: no cover
        return self.dict.keys()

    def values(self):  # pragma: no cover
        return self.dict.values()

    def items(self):  # pragma: no cover
        return self.dict.items()

    def iterkeys(self):  # pragma: no cover
        return self.dict.iterkeys()

    def itervalues(self):  # pragma: no cover
        for value in self.dict.itervalues():
            yield value

    def iteritems(self):  # pragma: no cover
        for key, value in self.dict.iteritems():
            yield key, value

    def update(self, *args, **kwds):  # pragma: no cover
        self.dict.update(*args, **kwds)

    def pop(self, key, defaukt=OrderedDict._OrderedDict__marker):  # pragma: no cover
        return self.dict.pop(key, default)

    def setdefault(self, key, default=None):  # pragma: no cover
        return self.dict.setdefault(key, default)

    def popitem(self, last=True):  # pragma: no cover
        return self.dict.popitem(last)

    def __reduce__(self):  # pragma: no cover
        return self.dict.__reduce__()

    def copy(self):  # pragma: no cover
        return self.dict.copy()

    def __eq__(self, other):  # pragma: no cover
        if isinstance(other, (self.__class__, OrderedDict)):
            return dict.__eq__(self, other) and all(_imap(_eq, self, other))
        return dict.__eq__(self, other)

    def __ne__(self, other):  # pragma: no cover
        return self.dict.__ne__(other)

    def viewkeys(self):  # pragma: no cover
        return self.dict.viewkeys()

    def viewvalues(self):  # pragma: no cover
        return self.dict.viewvalues()

    def viewitems(self):  # pragma: no cover
        return self.dict.viewitems()

    @classmethod
    def fromkeys(cls, iterable, value=None):  # pragma: no cover
        return self.dict.__class__.fromkeys(iterable, value)


class LRUCache(BaseCache):
    @property
    def to_be_removed(self):
        return self.iteritems().next()

    def remove(self):
        """Remove the least recently used item."""
        self.callback(*self.popitem(last=False))

    def __getitem__(self, key):
        print self.keys()
        if key in self:
            value = self.dict.__getitem__(key)
            del self[key]
            self.dict.__setitem__(key, value)
            return value
        raise KeyError

    def __setitem__(self, key, value):
        if key in self:
            del self[key]
        if self.is_full:
            self.remove()
        self.dict.__setitem__(key, value)


class LFUCache(BaseCache):
    def __init__(self, maxsize=float('inf')):
        self.usage_counter = Counter()
        BaseCache.__init__(self, maxsize)

    @property
    def to_be_removed(self):
        min_ = float('inf')
        lfu_key = None
        for k, v in self.usage_counter.iteritems():
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
        value = self.dict.__getitem__(key)
        self.usage_counter[key] += 1
        return value

    def __setitem__(self, key, value):
        self.usage_counter[key] += 1
        if self.is_full:
            self.remove()
        self.dict.__setitem__(key, value)
