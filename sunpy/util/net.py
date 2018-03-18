# -*- coding: utf-8 -*-
# Author: Florian Mayer <florian.mayer@bitsrc.org>

from __future__ import absolute_import, division, print_function
import os
import re
import sys
import shutil

# For Content-Disposition parsing
from sunpy.extern.six.moves.urllib.parse import urlparse, urljoin
from sunpy.extern.six.moves.urllib.request import urlopen
from sunpy.extern.six.moves.urllib.error import HTTPError, URLError
from sunpy.extern.six.moves import filter

from email.parser import FeedParser
from unicodedata import normalize

from sunpy.util import replacement_filename
from sunpy.extern import six

__all__ = ['slugify', 'get_content_disposition', 'get_filename',
           'get_system_filename', 'get_system_filename_slugify',
           'download_file', 'download_fileobj', 'check_download_file',
           'url_exists']

# Characters not allowed in slugified version.
_punct_re = re.compile(r'[:\t !"#$%&\'()*\-/<=>?@\[\\\]^_`{|},.]+')


def slugify(text, delim=u'_', encoding="ascii"):
    """ Slugify given unicode text. """
    text = normalize('NFKD', text)

    period = u'.'

    name_and_extension = text.rsplit(period, 1)
    name = name_and_extension[0]

    name = six.text_type(delim).join(
        filter(None, (word for word in _punct_re.split(name.lower()))))

    if len(name_and_extension) == 2:
        extension = name_and_extension[1]
        return six.text_type(period).join([name, extension])
    else:
        return name


def get_content_disposition(content_disposition):
    """ Get content disposition filename from given header. Do not include
    "Content-Disposition:". Returns a unicode string! """
    parser = FeedParser()
    parser.feed('Content-Disposition: ' + content_disposition)
    name = parser.close().get_filename()
    if not isinstance(name, six.text_type):
        name = name.decode('latin1', 'ignore')
    return name


def get_filename(sock, url):
    """ Get filename from given urllib2.urlopen object and URL.
    First, tries Content-Disposition, if unavailable, extracts
    name from URL. """
    name = None
    # NOTE: This gives bytes on 2 and unicode on 3.
    # How does 3.x know the encoding?
    cd = sock.headers.get('Content-Disposition', None)
    if cd is not None:
        try:
            name = get_content_disposition(cd)
        except IndexError:
            pass

    if not name:
        parsed = urlparse(url)
        name = parsed.path.rstrip('/').rsplit('/', 1)[-1]
    return six.text_type(name)


def get_system_filename(sock, url, default=u"file"):
    """ Get filename from given urllib2.urlopen object and URL.
    First, attempts to extract Content-Disposition, second, extract
    from URL, eventually fall back to default. Returns bytestring
    in file system encoding. """
    name = get_filename(sock, url)
    if not name:
        name = default
    return name.encode(sys.getfilesystemencoding(), 'ignore')


def get_system_filename_slugify(sock, url, default=u"file"):
    """ Get filename from given urllib2.urlopen object and URL.
    First, attempts to extract Content-Disposition, second, extract
    from URL, eventually fall back to default. Returns bytestring
    in file system encoding, normalized so it shouldn't violate
    operating system restrictions. """
    return slugify(get_system_filename(sock, url, default))


def download_file(url, directory, default=u'file', overwrite=False):
    """ Download file from url into directory. Try to get filename from
    Content-Disposition header, otherwise get from path of url. Fall
    back to default if both fail. Only overwrite existing files when
    overwrite is True. """
    opn = urlopen(url)
    try:
        path = download_fileobj(opn, directory, url, default, overwrite)
    finally:
        opn.close()
    return path


def download_file_with_cache(remote_url, directory, cache=True,
                             show_progress=True, timeout=None, overwrite=False,
                             default=u"file"):
    """
    Accepts a URL, downloads and optionally caches the result
    returning the filename, with a name determined by the file's MD5
    hash. If ``cache=True`` and the file is present in the cache, just
    returns the filename.

    Parameters
    ----------
    remote_url : str
        The URL of the file to download

    cache : bool, optional
        Whether to use the cache

    show_progress : bool, optional
        Whether to display a progress bar during the download (default
        is `True`)

    timeout : float, optional
        The timeout, in seconds.  Otherwise, use
        `astropy.utils.data.Conf.remote_timeout`.

    Returns
    -------
    local_path : str
        Returns the local path that the file was download to.

    Raises
    ------
    urllib2.URLError, urllib.error.URLError
        Whenever there's a problem getting the remote file.
    """
    from astropy.utils.console import ProgressBarOrSpinner
    from astropy.utils.data import (Conf, _get_download_cache_locs,
                                    _acquire_download_cache_lock,
                                    _release_download_cache_lock,
                                    check_free_space_in_dir)

    import urllib.error
    import shelve
    import socket

    conf = Conf()
    _dataurls_to_alias = {}

    if timeout is None:
        timeout = conf.remote_timeout

    missing_cache = False

    if cache:
        try:
            _, urlmapfn = _get_download_cache_locs()
        except OSError as e:
            msg = 'Remote data cache could not be accessed due to '
            estr = '' if len(e.args) < 1 else (': ' + str(e))
            warn(CacheMissingWarning(msg + e.__class__.__name__ + estr))
            cache = False
            # indicates that the cache is missing to raise a warning later
            missing_cache = True
    url_key = remote_url
    opn = urlopen(remote_url)
    filename = get_system_filename(opn, remote_url, default)
    path = os.path.join(directory, filename.decode('utf-8'))
    # Check if URL is Astropy data server, which has alias, and cache it.
    if (url_key.startswith(conf.dataurl) and
            conf.dataurl not in _dataurls_to_alias):
        with urllib.request.urlopen(conf.dataurl, timeout=timeout) as remote:
            _dataurls_to_alias[conf.dataurl] = [conf.dataurl, remote.geturl()]

    try:
        if cache:
            # We don't need to acquire the lock here, since we are only reading
            with shelve.open(urlmapfn) as url2hash:
                if url_key in url2hash:
                    if url2hash[url_key] != path:
                        shutil.move(url2hash[url_key], path)
                        url2hash[url_key] = path
                    return url2hash[url_key]
                # If there is a cached copy from mirror, use it.
                else:
                    for cur_url in _dataurls_to_alias.get(conf.dataurl, []):
                        if url_key.startswith(cur_url):
                            url_mirror = url_key.replace(cur_url,
                                                         conf.dataurl_mirror)
                            if url_mirror in url2hash:
                                return url2hash[url_mirror]
        with urllib.request.urlopen(remote_url, timeout=timeout) as remote:
            info = remote.info()
            if 'Content-Length' in info:
                try:
                    size = int(info['Content-Length'])
                except ValueError:
                    size = None
            else:
                size = None

            if size is not None:
                check_free_space_in_dir(directory, size)

            if show_progress:
                progress_stream = sys.stdout
            else:
                progress_stream = io.StringIO()

            dlmsg = "Downloading {0}".format(remote_url)
            with ProgressBarOrSpinner(size, dlmsg, file=progress_stream) as p:
                if not overwrite and os.path.exists(path):
                    path = replacement_filename(path)
                with open(path, 'wb') as f:
                    try:
                        bytes_read = 0
                        block = remote.read(conf.download_block_size)
                        while block:
                            f.write(block)
                            bytes_read += len(block)
                            p.update(bytes_read)
                            block = remote.read(conf.download_block_size)
                    except BaseException:
                        if os.path.exists(f.name):
                            os.remove(f.name)
                        raise

        if cache:
            _acquire_download_cache_lock()
            try:
                with shelve.open(urlmapfn) as url2hash:
                    # We check now to see if another process has
                    # inadvertently written the file underneath us
                    # already
                    if url_key in url2hash:
                        return url2hash[url_key]
                    url2hash[url_key] = f.name
            finally:
                _release_download_cache_lock()

    except urllib.error.URLError as e:
        if hasattr(e, 'reason') and hasattr(e.reason, 'errno') and e.reason.errno == 8:
            e.reason.strerror = e.reason.strerror + '. requested URL: ' + remote_url
            e.reason.args = (e.reason.errno, e.reason.strerror)
        raise e
    except socket.timeout as e:
        # this isn't supposed to happen, but occasionally a socket.timeout gets
        # through.  It's supposed to be caught in `urrlib2` and raised in this
        # way, but for some reason in mysterious circumstances it doesn't. So
        # we'll just re-raise it here instead
        raise urllib.error.URLError(e)
    return f.name


def download_fileobj(opn, directory, url='', default=u"file", overwrite=False):
    """ Download file from url into directory. Try to get filename from
    Content-Disposition header, otherwise get from path of url if given.
    Fall back to default if both fail. Only overwrite existing files when
    overwrite is True. """
    filename = get_system_filename(opn, url, default)
    path = os.path.join(directory, filename.decode('utf-8'))
    if not overwrite and os.path.exists(path):
        path = replacement_filename(path)
    with open(path, 'wb') as fd:
        shutil.copyfileobj(opn, fd)
    return path


def check_download_file(filename, remotepath, download_dir, remotename=None,
                        replace=False):
    """
    Downloads a file from remotepath to localpath if it isn't there.

    This function checks whether a file with name filename exists in the
    location, localpath, on the user's local machine.  If it doesn't,
    it downloads the file from remotepath.

    Parameters
    ----------
    filename : string
        Name of file.

    remotepath : string
        URL of the remote location from which filename can be downloaded.

    download_dir : string
        The files directory.

    remotename : (optional) string
        filename under which the file is stored remotely.
        Default is same as filename.

    replace : (optional) bool
        If True, file will be downloaded whether or not file already exists
        locally.

    Examples
    --------
    >>> from sunpy.util.net import check_download_file
    >>> remotepath = "http://www.download_repository.com/downloads/"
    >>> check_download_file("filename.txt", remotepath, download_dir='.')   # doctest: +SKIP
    """
    # Check if file already exists locally.  If not, try downloading it.
    if replace or not os.path.isfile(os.path.join(download_dir, filename)):
        # set local and remote file names be the same unless specified
        # by user.
        if not isinstance(remotename, six.string_types):
            remotename = filename

        download_file(urljoin(remotepath, remotename),
                      download_dir, default=filename, overwrite=replace)


def url_exists(url, timeout=2):
    """
    Checks whether a url is online.

    Parameters
    ----------
    url: `str`
        A string containing a URL

    Returns
    -------
    value: `bool`

    Examples
    --------
    >>> from sunpy.util.net import url_exists
    >>> url_exists('http://www.google.com')  #doctest: +REMOTE_DATA
    True
    >>> url_exists('http://aslkfjasdlfkjwerf.com')  #doctest: +REMOTE_DATA
    False
    """
    try:
        urlopen(url, timeout=timeout)
    except HTTPError:
        return False
    except URLError:
        return False
    else:
        return True


def is_online():
    """
    Checks whether an internet connection is available.

    Returns
    -------
    value: `bool`

    Examples
    --------
    >>> from sunpy.util.net import is_online
    >>> is_online()  #doctest: +REMOTE_DATA
    True
    """
    return url_exists('http://www.google.com')
