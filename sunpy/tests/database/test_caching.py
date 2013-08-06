from __future__ import absolute_import

from collections import deque

import pytest

from sunpy.database.caching import BaseCache, LRUCache, LFUCache


def test_custom_cache():
    class FIFO(BaseCache):
        def __init__(self, maxsize=float('inf')):
            self.queue = deque([], maxsize)
            BaseCache.__init__(self, maxsize)

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
    lrucache[1]
    lrucache[3]
    lrucache[4] = 'd'
    assert len(lrucache) == 3
    assert lrucache[1] == 'a'
    assert lrucache[3] == 'c'
    assert lrucache[4] == 'd'
    #assert lrucache.items() == [(3, 'c'), (1, 'a'), (4, 'd')]


def test_lru_cache_missing_item():
    with pytest.raises(KeyError):
        LRUCache()[0]


def test_lfu_cache():
    lfucache = LFUCache(3)
    lfucache[1] = 'a'
    lfucache[2] = 'b'
    lfucache[3] = 'c'
    lfucache[1]
    lfucache[2]
    lfucache[4] = 'd'
    assert len(lfucache) == 3
    assert lfucache[1] == 'a'
    assert lfucache[2] == 'b'
    assert lfucache[4] == 'd'
    #assert lrucache.items() == [(1, 'a'), (2, 'b'), (4, 'd')]
