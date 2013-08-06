from __future__ import absolute_import

from abc import ABCMeta, abstractmethod
from collections import OrderedDict, Counter

__all__ = ['BaseCache', 'LRUCache', 'LFUCache']


class BaseCache(OrderedDict):
    """BaseCache is a class that inherits from OrderedDict. It has a
    certain capacity, stored in the attribute :attr:`maxsize`. Whether this
    capacity is reached, can be checked by using the boolean property
    :attr:`is_full`. To implement a custom cache, inherit from this class and
    override the methods ``__getitem__`` and ``__setitem__``.
    Call the method :meth:`callback` as soon as an item from the cache is
    removed.

    Methods
    -------
    callback

    """
    __metaclass__ = ABCMeta

    def __init__(self, maxsize=float('inf')):
        self.maxsize = maxsize
        OrderedDict.__init__(self)

    def set_size(self, size):
        self.maxsize = size
        while len(self) > self.maxsize:
            self.remove()

    @abstractmethod
    def __getitem__(self, key):
        return  # pragma: no cover

    @abstractmethod
    def __setitem__(self, key, value):
        return  # pragma: no cover

    @abstractmethod
    def remove(self):
        """Call this method to manually remove one item from the cache. Which
        item is removed, depends on the implementation of the cache. After the
        item has been removed, the callback method is called.

        """
        return  # pragma: no cover

    def callback(self, key, value):
        """This method should be called (by convention) if an item is removed
        from the cache because it is full. The passed key and value are the
        ones that are removed. By default this method does nothing, but it
        can be customized in a custom cache that inherits from this base class.

        """

    @property
    def is_full(self):
        return len(self) == self.maxsize


class LRUCache(BaseCache):
    def remove(self):
        """Remove the least recently used item."""
        self.callback(*self.popitem(last=False))

    def __getitem__(self, key):
        if key in self:
            value = OrderedDict.__getitem__(self, key)
            del self[key]
            OrderedDict.__setitem__(self, key, value)
            return value
        raise KeyError

    def __setitem__(self, key, value):
        if key in self:
            del self[key]
        if self.is_full:
            self.remove()
        OrderedDict.__setitem__(self, key, value)


class LFUCache(BaseCache):
    def __init__(self, maxsize=float('inf')):
        self.usage_counter = Counter()
        BaseCache.__init__(self, maxsize)

    def remove(self):
        """Remove the least frequently used item."""
        lfu_key = self._get_lfu_key()
        val = self[lfu_key]
        del self[lfu_key]
        del self.usage_counter[lfu_key]
        self.callback(lfu_key, val)

    def _get_lfu_key(self):
        min_ = float('inf')
        for k, v in self.usage_counter.iteritems():
            if v < min_:
                min_ = v
                key = k
        return key

    def __getitem__(self, key):
        value = OrderedDict.__getitem__(self, key)
        self.usage_counter[key] += 1
        return value

    def __setitem__(self, key, value):
        self.usage_counter[key] += 1
        if self.is_full:
            self.remove()
        OrderedDict.__setitem__(self, key, value)
