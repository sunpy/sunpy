# -*- coding: utf-8 -*-
# Author: Florian Mayer <florian.mayer@bitsrc.org>

#pylint: disable=W0613

from __future__ import absolute_import

import pytest

import os
import tempfile

from functools import partial

from sunpy.net.download import Downloader, default_name
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
    th = threading.Thread(target=dw.wait)
    th.daemon = True
    th.start()
    x.wait(10)
    assert x.isSet()
    dw.stop()


def test_download_files():
    items = []

    def wait_for(n, callback):  # pylint: disable=W0613
        def _fun(handler):
            items.append(handler)
            if len(items) == n:
                callback(items)
        return _fun

    tmp = tempfile.mkdtemp()
    path_fun = partial(default_name, tmp)

    dw = Downloader(1, 1)

    on_finish = wait_for(2, lambda _: dw.stop())
    dw.download('ftp://speedtest.inode.at/speedtest-5mb', path_fun, on_finish)
    dw.download('ftp://speedtest.inode.at/speedtest-20mb', path_fun, on_finish)
    # dw.download('ftp://speedtest.inode.at/speedtest-100mb', path_fun, on_finish)

    dw.wait()

    for item in items:
        assert os.path.exists(item ['path'])
