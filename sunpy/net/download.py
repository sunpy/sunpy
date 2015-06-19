# -*- coding: utf-8 -*-
# Author: Florian Mayer <florian.mayer@bitsrc.org>
#
# This module was developed with funding provided by
# the ESA Summer of Code (2011).


from __future__ import absolute_import

import os
import re
import urllib2
import threading

from functools import partial
from contextlib import closing
from collections import defaultdict, deque

import sunpy as spy

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

            with closing(urllib2.urlopen(url)) as sock:
                fullname = path(sock, url)

                with open(fullname, 'wb') as fd:
                    while True:
                        rec = sock.read(self.buf)
                        if not rec:
                            with self.mutex:
                                self._close(callback, [{'path': fullname}], server)
                            break
                        else:
                            fd.write(rec)
        except Exception, e:
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
        self.done_lock.acquire()

    def stop(self):
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
        # @todo: explain

        server = self._get_server(url)

        # Create function to compute the filepath to download to if not set
        default_dir = spy.config.get("downloads", "download_dir")

        if path is None:
            path = partial(default_name, default_dir)
        elif isinstance(path, basestring):
            path = partial(default_name, path)

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
            for k, v in self.q.iteritems():  # pylint: disable=W0612
                while v:
                    if self._attempt_download(*v[0]):
                        v.popleft()
                        if self.conns == self.max_total:
                            return
                    else:
                        break
