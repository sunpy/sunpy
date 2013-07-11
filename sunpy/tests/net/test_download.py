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


def test_download_http():
    items = []
    lck = threading.Lock()

    def wait_for(n, callback):  # pylint: disable=W0613
        def _fun(handler):
            with lck:
                items.append(handler)
                if len(items) == n:
                    callback(items)
        return _fun

    tmp = tempfile.mkdtemp()
    path_fun = partial(default_name, tmp)

    dw = Downloader(1, 1)
    # If this fires, the assertion below will fail.
    threading.Timer(60, dw.stop)

    on_finish = wait_for(3, lambda _: dw.stop())
    dw.download('http://ajax.googleapis.com/ajax/libs/jquery/2.0.0/jquery.min.js', path_fun, on_finish)
    dw.download('http://ajax.googleapis.com/ajax/libs/webfont/1.4.2/webfont.js', path_fun, on_finish)
    dw.download('https://raw.github.com/sunpy/sunpy/master/INSTALL.txt', path_fun, on_finish)
    # dw.download('ftp://speedtest.inode.at/speedtest-100mb', path_fun, on_finish)

    dw.wait()

    assert len(items) == 3

    for item in items:
        assert os.path.exists(item ['path'])
