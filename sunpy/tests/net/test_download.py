# -*- coding: utf-8 -*-
# Author: Florian Mayer <florian.mayer@bitsrc.org>

#pylint: disable=W0613

from __future__ import absolute_import

import pytest

import os
import tempfile
import threading

from functools import partial

import sunpy

from sunpy.net.download import Downloader, default_name


class CalledProxy(object):
    def __init__(self, fn):
        self.fn = fn
        self.fired = False

    def __call__(self, *args, **kwargs):
        self.fn(*args, **kwargs)
        self.fired = True


class MockConfig(object):
    def __init__(self):
        self.dct = {}

    def add_section(self, name, dct):
        self.dct[name] = dct

    def get(self, one, other):
        return self.dct[one][other]


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
    _stop = lambda _: dw.stop()

    timeout = CalledProxy(_stop)
    threading.Timer(60, timeout).start()

    on_finish = wait_for(3, lambda _: dw.stop())
    dw.download('http://ajax.googleapis.com/ajax/libs/jquery/2.0.0/jquery.min.js', path_fun, on_finish)
    dw.download('http://ajax.googleapis.com/ajax/libs/webfont/1.4.2/webfont.js', path_fun, on_finish)
    dw.download('https://raw.github.com/sunpy/sunpy/master/INSTALL.txt', path_fun, on_finish)
    # dw.download('ftp://speedtest.inode.at/speedtest-100mb', path_fun, on_finish)

    dw.wait()

    assert len(items) == 3
    assert not timeout.fired

    for item in items:
        assert os.path.exists(item['path'])


def test_download_default_dir():
    _config = sunpy.config

    try:
        tmpdir = tempfile.mkdtemp()

        sunpy.config = MockConfig()
        sunpy.config.add_section(
            "downloads", {"download_dir": tmpdir}
        )

        dw = Downloader(1, 1)
        _stop = lambda _: dw.stop()

        timeout = CalledProxy(_stop)
        errback = CalledProxy(_stop)
        dw.download(
            'http://ajax.googleapis.com/ajax/libs/jquery/2.0.0/jquery.min.js',
            callback=_stop,
            errback=errback
        )

        threading.Timer(10, timeout).start()
        dw.wait()

        assert not timeout.fired
        assert not errback.fired
        assert os.path.exists(os.path.join(tmpdir, 'jquery.min.js'))
    finally:
        sunpy.config = _config

def test_download_dir():
    tmpdir = tempfile.mkdtemp()

    dw = Downloader(1, 1)
    _stop = lambda _: dw.stop()
    timeout = CalledProxy(_stop)
    errback = CalledProxy(_stop)

    dw.download(
        'http://ajax.googleapis.com/ajax/libs/jquery/2.0.0/jquery.min.js',
        tmpdir,
        callback=_stop,
        errback=errback
    )

    threading.Timer(10, timeout).start()
    dw.wait()
    assert not timeout.fired
    assert not errback.fired
    assert os.path.exists(os.path.join(tmpdir, 'jquery.min.js'))
