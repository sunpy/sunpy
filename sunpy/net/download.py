# -*- coding: utf-8 -*-
# Author: Florian Mayer <florian.mayer@bitsrc.org>
#
# This module was developed with funding provided by
# the ESA Summer of Code (2011).


from __future__ import absolute_import

import os
import re
import urllib2
import select
import socket
import threading

from functools import partial
from collections import defaultdict, deque

import sunpy

class IDPool(object):
    """
    Pool that returns unique identifiers in a thread-safe way.
    
    Identifierers obtained using the get method are guaranteed to not be
    returned by it again until they are released using the release method.
    
        >>> pool = IDPool()
        >>> pool.get()
        0
        >>> pool.get()
        1
        >>> pool.get()
        2
        >>> pool.release(1)
        >>> pool.get()
        1
        >>> 
    """
    def __init__(self):
        self.max_id = -1
        self.free_ids = []
        
        self._lock = threading.Lock()
    
    def get(self):
        """ Return a new integer that is unique in this pool until
        it is released. """
        self._lock.acquire()
        try:
            if self.free_ids:
                return self.free_ids.pop()
            else:
                self.max_id += 1
                return self.max_id
        finally:
            self._lock.release()
    
    def release(self, id_):
        """ Release the id. It can now be returned by get again.
        
        Will reset the IDPool if the last id in use is released. """
        self._lock.acquire()
        try:
            self.free_ids.append(id_)
            if len(self.free_ids) == self.max_id + 1:
                self.reset()
        finally:
            self._lock.release()
    
    def reset(self):
        """ Reset the state of the IDPool. This should only be called when
        no identifier is in use. """
        self.max_id = -1
        self.free_ids = []


def socketpair():
    """ Return pair of connected sockets. Unlike socket.socketpair this
    is platform independant. However, if socket.socketpair is available,
    it is used here as well. """
    if hasattr(socket, 'socketpair'):
        # Unix.
        return socket.socketpair()
    
    try:
        acceptor = socket.socket()
        # Random port. Only accept local connections.
        acceptor.bind(('127.0.0.1', 0))
        # We know we'll only get one connection.
        acceptor.listen(1)

        one = socket.socket()
        one.connect(acceptor.getsockname())
        
        other = acceptor.accept()[0]
    finally:
        acceptor.close()
    return one, other


class Reactor(object):
    def __init__(self):
        self.syncr, self.synce = socketpair()
        self.ids = IDPool()
        self.call = []
        self.tcalls = {}
        self.callb = {}
        self.running = False
    
    def _unset_running(self):
        self.running = False
    
    def stop(self):
        self.call_sync(self._unset_running)
    
    def run(self, tcallfreq=0):
        self.running = True
        while self.running:
            ret = self.poll(tcallfreq if self.tcalls else None)
            self._call_calls()
            self._call_tcalls()
            for fd in ret:
                try:
                    fun = self.callb[fd]
                except KeyError:
                    continue
                fun()
    
    def poll(self, timeout=None):
        raise NotImplementedError
    
    def call_sync(self, fun):
        self.call.append(fun)
        self.synce.send('m')
    
    def _call_tcalls(self):
        for fun in list(self.tcalls.itervalues()):
            fun()
    
    def get_tid(self):
        return self.ids.get()
    
    def add_tcall(self, id_, fun):
        self.tcalls[id_] = fun
        return id_
    
    def rem_tcall(self, id_):
        del self.tcalls[id_]
        self.ids.release(id_)
    
    def _call_calls(self):
        for fun in self.call:
            fun()
            self.syncr.recv(1)
        self.call = []
    
    def add_fd(self, fd, callback):
        self.callb[fd] = callback
    
    def rem_fd(self, fd):
        del self.callb[fd]


class DummyReactor(Reactor): 
    def run(self, tcallfreq=0):
        self.running = True
        while self.running:
            select.select(
                [self.syncr], [], [], tcallfreq if self.tcalls else None
            )
            self._call_calls()
            self._call_tcalls()
    

class SelectReactor(Reactor):
    avail = hasattr(select, 'select')
    def __init__(self, fds=None):
        Reactor.__init__(self)
        self.fds = set() if fds is None else fds
        self.fds.add(self.syncr)
    
    def poll(self, timeout=None):
        return select.select(self.fds, [], [], timeout)[0]
    
    def add_fd(self, fd, callback):
        super(SelectReactor, self).add_fd(fd, callback)
        self.fds.add(fd)
        
    def rem_fd(self, fd):
        super(SelectReactor, self).rem_fd(fd)
        self.fds.remove(fd)
    
    def close(self):
        pass


class PollReactor(Reactor):
    avail = hasattr(select, 'poll')
    def __init__(self, fds=None):
        Reactor.__init__(self)
        self.poller = select.poll()
        
        self.poller.register(
            self.syncr,
            select.POLLERR | select.POLLHUP | select.POLLNVAL | select.POLLIN
        )
    
    def poll(self, timeout=None):
        if timeout is not None:
            timeout = int(timeout * 1000)
        return (fileno for fileno, flags in self.poller.poll(timeout))
    
    def add_fd(self, fd, callback):
        fd = fd.fileno()
        
        super(PollReactor, self).add_fd(fd, callback)
        self.poller.register(
            fd,
            select.POLLERR | select.POLLHUP | select.POLLNVAL | select.POLLIN
        )
        
    def rem_fd(self, fd):
        fd = fd.fileno()
        
        super(PollReactor, self).rem_fd(fd)
        self.poller.unregister(fd)
    
    def close(self):
        self.poller.close()


class EPollReactor(Reactor):
    avail = hasattr(select, 'epoll')
    def __init__(self, fds=None):
        Reactor.__init__(self)
        self.poller = select.epoll()
        
        self.poller.register(
            self.syncr,
            select.POLLERR | select.POLLHUP | select.POLLNVAL | select.POLLIN
        )
    
    def poll(self, timeout=None):
        if timeout is None:
            timeout = -1
        return (fileno for fileno, flags in self.poller.poll(timeout))
    
    def add_fd(self, fd, callback):
        fd = fd.fileno()
        
        super(EPollReactor, self).add_fd(fd, callback)
        self.poller.register(
            fd,
            select.EPOLLERR | select.EPOLLHUP | select.POLLNVAL | select.EPOLLIN
        )
        
    def rem_fd(self, fd):
        fd = fd.fileno()
        
        super(EPollReactor, self).rem_fd(fd)
        self.poller.unregister(fd)
    
    def close(self):
        self.poller.close()


# FIXME: Add KQueueReactor.
DReactor = None
for reactor in [SelectReactor, PollReactor, EPollReactor]:
    if reactor.avail:
        DReactor = reactor

if DReactor is None:
    # This really should not be happening.
    raise EnvironmentError('No suitable function in select module.')


def default_name(path, sock, url):
    name = sock.headers.get('Content-Disposition', url.rsplit('/', 1)[-1])
    return os.path.join(path, name)


class Downloader(object):
    def __init__(self, max_conn=5, max_total=20):
        self.max_conn = max_conn
        self.max_total = max_total
        self.conns = 0
        
        self.connections = defaultdict(int) # int() -> 0
        self.q = defaultdict(deque)
        
        self.reactor = DummyReactor()
        self.buf = 9096
    
    def _download(self, sock, fd, callback, id_=None):
        rec = sock.read(self.buf)
        if not rec:
            callback()
            
            if id_ is not None:
                self.reactor.rem_tcall(id_)
            else:
                self.reactor.rem_fd(sock)
            fd.close()
        else:
            fd.write(rec)
    
    def _start_download(self, url, path, callback, errback):
        try:
            server = self._get_server(url)
            
            self.connections[server] += 1
            self.conns += 1
            
            sock = urllib2.urlopen(url)
            fullname = path(sock, url)
            
            args = [
                sock, open(fullname, 'wb'),
                partial(self._close, callback, [{'path': fullname}], server),
            ]
            
            try:
                # hasattr does not work because HTTPResponse objects have a
                # fileno method that raises AttributeError when called.
                # Don't ask me.
                sock.fileno()
                # Due to caching it's not save to use select.
                raise AttributeError
            except AttributeError:
                id_ = self.reactor.get_tid()
                self.reactor.add_tcall(
                    id_, partial(self._download, *args + [id_])
                )
            else:
                self.reactor.add_fd(sock, partial(self._download, *args))
        except Exception, e:
            if errback is not None:
                errback(e)
    
    def _attempt_download(self, url, path, callback, errback):
        num_connections = self.connections[self._get_server(url)]
        
        # If max downloads has not been exceeded, begin downloading
        if (num_connections < self.max_conn and self.conns < self.max_total):
            self._start_download(url, path, callback, errback)
            return True
        return False

    def _get_server(self, url):
        """Returns the server name for a given URL.
        
        Examples: http://server.com, server.org, ftp.server.org, etc.
        """
        return re.search('(\w+://)?([\w\.]+)', url).group(2)
        
    def _default_callback(self, *args):
        """Default callback to execute on a successfull download"""
        pass
        
    def _default_error_callback(self, e):
        """Default callback to execute on a failed download"""
        raise e

    def download(self, url, path=None, callback=None, errback=None):
        self.reactor.call_sync(
            lambda: self._dwn(url, path, callback, errback)
        )

    def _dwn(self, url, path=None, callback=None, errback=None):
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
        default_dir = sunpy.config.get("downloads", "download_dir")

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
        callback(*args)
        
        if self.q[server]:
            self._start_download(*self.q[server].pop())
        else:
            self.connections[server] -= 1
            self.conns -= 1
            
            for k, v in self.q.iteritems(): #pylint: disable=W0612
                while v:
                    if self._attempt_download(*v[0]):
                        v.popleft()
                        if self.conns == self.max_total:
                            return
                    else:
                        break


if __name__ == '__main__':
    import tempfile

    def wait_for(n, callback): #pylint: disable=W0613
        items = []
        def _fun(handler):
            items.append(handler)
            if len(items) == n:
                callback(items)
        return _fun
    
    
    tmp = tempfile.mkdtemp()
    print tmp
    path_fun = partial(default_name, tmp)
    
    dw = Downloader(1, 2)
    
    on_finish = wait_for(4, lambda _: dw.reactor.stop())
    dw.download('ftp://speedtest.inode.at/speedtest-5mb', path_fun, on_finish)
    dw.download('ftp://speedtest.inode.at/speedtest-20mb', path_fun, on_finish)
    dw.download('https://bitsrc.org', path_fun, on_finish)
    dw.download('ftp://speedtest.inode.at/speedtest-100mb', path_fun, on_finish)
    
    print dw.conns
    
    dw.reactor.run()
