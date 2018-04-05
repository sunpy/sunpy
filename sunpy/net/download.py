# -*- coding: utf-8 -*-
# Author: Florian Mayer <florian.mayer@bitsrc.org>
#
# This module was developed with funding provided by
# the ESA Summer of Code (2011).


from __future__ import absolute_import

import os
import re
import threading

from functools import partial
from contextlib import closing
from collections import defaultdict, deque

from sunpy.extern import six
from sunpy.extern.six.moves import urllib
from sunpy.extern.six import iteritems

import sunpy
from sunpy.util.progressbar import TTYProgressBar as ProgressBar
from sunpy.util.config import get_and_create_download_dir


__all__  = ['Downloader', 'Results']

def default_name(path, sock, url):
    name = sock.headers.get('Content-Disposition', url.rsplit('/', 1)[-1])
    return os.path.join(path, name)


class Downloader(object):
    def __init__(self, max_conn=5, max_total=20):
        self.max_conn = max_conn
        self.max_total = max_total
        self.conns = 0

        self.connections = defaultdict(int)  # int() -> 0
        self.q = defaultdict(deque)

        self.buf = 9096

        self.done_lock = threading.Semaphore(0)
        self.mutex = threading.Lock()

    def _start_download(self, url, path, callback, errback):
        try:
            server = self._get_server(url)

            with self.mutex:
                self.connections[server] += 1
                self.conns += 1

            with closing(urllib.request.urlopen(url)) as sock:
                fullname = path(sock, url)
                dir_ = os.path.abspath(os.path.dirname(fullname))
                if not os.path.exists(dir_):
                    os.makedirs(dir_)

                with open(fullname, 'wb') as fd:
                    while True:
                        rec = sock.read(self.buf)
                        if not rec:
                            with self.mutex:
                                self._close(callback, [{'path': fullname}], server)
                            break
                        else:
                            fd.write(rec)
        except Exception as e:
            # TODO: Fix the silent failing
            if errback is not None:
                with self.mutex:
                    self._close(errback, [e], server)

    def _attempt_download(self, url, path, callback, errback):
        """ Attempt download. If max. connection limit reached, queue for download later.
        """

        num_connections = self.connections[self._get_server(url)]

        # If max downloads has not been exceeded, begin downloading
        if num_connections < self.max_conn and self.conns < self.max_total:
            th = threading.Thread(
                target=partial(self._start_download, url,
                               path, callback, errback)
            )
            th.daemon = True
            th.start()
            return True
        return False

    def _get_server(self, url):
        """Returns the server name for a given URL.

        Examples: http://server.com, server.org, ftp.server.org, etc.
        """
        return re.search('(\w+://)?([\w\.]+)', url).group(2)

    def _default_callback(self, *args):
        """Default callback to execute on a successful download"""
        pass

    def _default_error_callback(self, e):
        """Default callback to execute on a failed download"""
        raise e

    def wait(self):
        """
        Waits for all files to download and then return.
        """
        self.done_lock.acquire()

    def stop(self):
        """
        Stops all downloads and then return.
        """
        self.done_lock.release()

    def init(self):
        pass

    def download(self, url, path=None, callback=None, errback=None):
        """Downloads a file at a specified URL.

        Parameters
        ----------
        url : string
            URL of file to download
        path : function, string
            Location to save file to. Can specify either a directory as a string
            or a function with signature: (path, url).
            Defaults to directory specified in sunpy configuration
        callback : function
            Function to call when download is successfully completed
        errback : function
            Function to call when download fails

        Returns
        -------
        out : None
        """
        # Load balancing?
        # TODO: explain

        server = self._get_server(url)

        # Create function to compute the filepath to download to if not set

        if path is None:
            path = partial(default_name, get_and_create_download_dir())
        elif isinstance(path, six.string_types):
            path = partial(default_name, path)
        elif not callable(path):
            raise ValueError("path must be: None, string or callable")

        # Use default callbacks if none were specified
        if callback is None:
            callback = self._default_callback
        if errback is None:
            errback = self._default_error_callback

        # Attempt to download file from URL
        if not self._attempt_download(url, path, callback, errback):
            # If there are too many concurrent downloads, queue for later
            self.q[server].append((url, path, callback, errback))

    def _close(self, callback, args, server):
        """ Called after download is done. Activated queued downloads, call callback.
        """
        callback(*args)

        self.connections[server] -= 1
        self.conns -= 1

        if self.q[server]:
            self._attempt_download(*self.q[server].pop())
        else:
            for k, v in iteritems(self.q):  # pylint: disable=W0612
                while v:
                    if self._attempt_download(*v[0]):
                        v.popleft()
                        if self.conns == self.max_total:
                            return
                    else:
                        break


class Results(object):
    """ Returned by VSOClient.get. Use .wait to wait
    for completion of download.
    """
    def __init__(self, callback, n=0, done=None):
        self.callback = callback
        self.n = self.total = n
        self.map_ = {}
        self.done = done
        self.evt = threading.Event()
        self.errors = []
        self.lock = threading.RLock()

        self.progress = None

    def submit(self, keys, value):
        """
        Submit

        Parameters
        ----------
        keys : list
            names under which to save the value
        value : object
            value to save
        """
        for key in keys:
            self.map_[key] = value
        self.poke()

    def poke(self):
        """ Signal completion of one item that was waited for. This can be
        because it was submitted, because it lead to an error or for any
        other reason. """
        with self.lock:
            self.n -= 1
            if self.progress is not None:
                self.progress.poke()
            if not self.n:
                if self.done is not None:
                    self.map_ = self.done(self.map_)
                self.callback(self.map_)
                self.evt.set()

    def require(self, keys):
        """ Require that keys be submitted before the Results object is
        finished (i.e., wait returns). Returns a callback method that can
        be used to submit the result by simply calling it with the result.

        keys : list
            name of keys under which to save the result
        """
        with self.lock:
            self.n += 1
            self.total += 1
            return partial(self.submit, keys)

    def wait(self, timeout=100, progress=True):
        """ Wait for result to be complete and return it. """
        # Giving wait a timeout somehow circumvents a CPython bug that the
        # call gets ininterruptible.
        if progress:
            with self.lock:
                self.progress = ProgressBar(self.total, self.total - self.n)
                self.progress.start()
                self.progress.draw()

        while not self.evt.wait(timeout):
            pass
        if progress:
            self.progress.finish()

        return self.map_

    def add_error(self, exception):
        """ Signal a required result cannot be submitted because of an
        error. """
        self.errors.append(exception)
        self.poke()
