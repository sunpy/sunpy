# -*- coding: utf-8 -*-
# Author: Florian Mayer <florian.mayer@bitsrc.org>

#pylint: disable=W0613

from __future__ import absolute_import

import pytest

from sunpy.net.download import Downloader
import threading

def wait_for(n, callback): #pylint: disable=W0613
    items = []
    def _fun(handler):
        items.append(handler)
        if len(items) == n:
            callback(items)
    return _fun


def path_fun(*args, **kwargs):
    raise ValueError


def test_path_exception():
    x = threading.Event()
    dw = Downloader(1, 2)
    dw.download(
        "http://google.at", path_fun, errback=wait_for(1, lambda a: x.set())
    )
    th = threading.Thread(target=dw.reactor.run)
    th.daemon = True
    th.start()
    assert x.wait(10)
    dw.reactor.stop()
