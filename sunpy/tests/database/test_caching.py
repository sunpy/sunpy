# Author: Simon Liedtke <liedtke.simon@googlemail.com>
#
# This module was developed with funding provided by
# the Google Summer of Code (2013).

from __future__ import absolute_import

from collections import deque

import pytest

from sunpy.database.caching import BaseCache, LRUCache, LFUCache


def test_custom_cache():
    class FIFO(BaseCache):
        def __init__(self, maxsize=float('inf')):
            self.queue = deque([], maxsize)
            BaseCache.__init__(self, maxsize)

        @property
        def to_be_removed(self):
            try:
                return self.queue[0]
            except IndexError:
                return None

        def remove(self):
            del self.queue[0]

        def __getitem__(self, key):
            for k, value in self.queue:
                if k == key:
                    return value
            raise KeyError

        def __setitem__(self, key, value):
            self.queue.append((key, value))

        def __len__(self):
            return len(self.queue)
    cache = FIFO(3)
    cache[1] = 'a'
    cache[2] = 'b'
    cache[3] = 'c'
    assert cache.to_be_removed == (1, 'a')
    cache[4] = 'd'
    assert len(cache) == 3
    assert cache[2] == 'b'
    assert cache[3] == 'c'
    assert cache[4] == 'd'


def test_lru_cache():
    lrucache = LRUCache(3)
    lrucache[1] = 'a'
    lrucache[2] = 'b'
    lrucache[3] = 'c'
    assert lrucache.to_be_removed == (1, 'a')
    lrucache[1]
    lrucache[3]
    assert lrucache.to_be_removed == (2, 'b')
    lrucache[4] = 'd'
    assert len(lrucache) == 3
    assert lrucache.to_be_removed == (1, 'a')
    assert lrucache == {1: 'a', 3: 'c', 4: 'd'}


def test_lru_cache_missing_item():
    with pytest.raises(KeyError):
        LRUCache()[0]


def test_lfu_cache():
    lfucache = LFUCache(3)
    lfucache[1] = 'a'
    lfucache[2] = 'b'
    lfucache[3] = 'c'
    assert lfucache.to_be_removed == (1, 'a')
    lfucache[1]
    lfucache[2]
    assert lfucache.to_be_removed == (3, 'c')
    lfucache[4] = 'd'
    assert len(lfucache) == 3
    assert lfucache.to_be_removed == (4, 'd')
    assert lfucache == {1: 'a', 2: 'b', 4: 'd'}
