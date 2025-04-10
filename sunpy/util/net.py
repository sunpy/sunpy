"""
This module provides general net utility functions.
"""
import os
import sys
import shutil
from unicodedata import normalize
from email.parser import FeedParser
from urllib.parse import urlparse
from urllib.request import urlopen

from sunpy.util import replacement_filename
from sunpy.util.decorators import deprecated

__all__ = ['parse_header', 'slugify', 'get_content_disposition', 'get_filename',
           'get_system_filename', 'download_file', 'download_fileobj']


def slugify(text, delim='_'):
    r"""
    Slugify given unicode text.

    This function performs a Unicode normalization to NFKC form, followed by replacing
    the following characters by the delimiter:

    : (tab) (space) ! " # $ % & ' ( ) * - / < = > ? @ [ \\ ] ^ _ ` { | } ,

    Parameters
    ----------
    text : `str`
        A `str` to slugify.
    delim : `str`, optional
        The delimiter for the input ``text``. Default is "_".

    Returns
    -------
    `str` :
        The slugify `str` name.
    """
    text = normalize('NFKC', text)

    chars_to_replace = ":\t !\"#$%&'()*-/<=>?@[\\]^_`{|},"
    trans_map = str.maketrans({c: delim for c in chars_to_replace})

    return text.translate(trans_map)


def get_content_disposition(content_disposition):
    """
    Get the content disposition filename from given header.

    **Do not include "Content-Disposition:".**

    Parameters
    ----------
    content_disposition : `str`
        The content disposition header.

    Returns
    -------
    `str` :
        The content disposition filename.
    """
    parser = FeedParser()
    parser.feed('Content-Disposition: ' + content_disposition)
    name = parser.close().get_filename()
    if name and not isinstance(name, str):
        name = name.decode('latin1', 'ignore')
    return name


def get_filename(sock, url):
    """
    Get filename from given `~urllib.request.urlopen` object and URL.

    First, tries the "Content-Disposition", if unavailable, extracts name from the URL.

    Parameters
    ----------
    sock : `~urllib.request.urlopen`
        The `~urllib.request.urlopen` to parse for the filename.
    url : `str`
        The URL to parse for the filename.

    Returns
    -------
    `str`
        The filename.
    """
    name = None
    cd = sock.headers.get('Content-Disposition', None)
    if cd is not None:
        try:
            name = get_content_disposition(cd)
        except IndexError:
            pass

    if not name:
        parsed = urlparse(url)
        name = parsed.path.rstrip('/').rsplit('/', 1)[-1]
    return str(name)


def get_system_filename(sock, url, default="file"):
    """
    Get filename from given `~urllib.request.urlopen` object and URL.

    First, tries the "Content-Disposition", if unavailable, extracts name from the URL.
    If this fails, the ``default`` keyword will be used.

    Parameters
    ----------
    sock : `~urllib.request.urlopen`
        The `~urllib.request.urlopen` to parse for the filename.
    url : `str`
        The URL to parse for the filename.
    default : `str`, optional
        The name to use if the first two methods fail. Defaults to "file".

    Returns
    -------
    `bytes`:
        The filename in file system encoding.
    """
    name = get_filename(sock, url)
    if not name:
        name = str(default)
    return name.encode(sys.getfilesystemencoding(), 'ignore')


@deprecated(since="7.0", alternative="parfive.Downloader.simple_download")
def download_fileobj(opn, directory, url='', default="file", overwrite=False):
    """
    Download a file from a url into a directory.

    Tries the "Content-Disposition", if unavailable, extracts name from the URL.
    If this fails, the ``default`` keyword will be used.

    Parameters
    ----------
    opn : `~urllib.request.urlopen`
        The `~urllib.request.urlopen` to download.
    directory : `str`
        The directory path to download the file in to.
    url : `str`
        The URL to parse for the filename.
    default : `str`, optional
        The name to use if the first two methods fail. Defaults to "file".
    overwrite : `bool`, optional
        If `True` will overwrite a file of the same name. Defaults to `False`.

    Returns
    -------
    `str`:
        The file path for the downloaded file.
    """
    filename = get_system_filename(opn, url, default)
    path = os.path.join(directory, filename.decode('utf-8'))
    if overwrite and os.path.exists(path):
        path = replacement_filename(path)
    with open(path, 'wb') as fd:
        shutil.copyfileobj(opn, fd)
    return path


@deprecated(since="7.0", alternative="parfive.Downloader.simple_download")
def download_file(url, directory, default="file", overwrite=False):
    """
    Download a file from a url into a directory.

    Tries the "Content-Disposition", if unavailable, extracts name from the URL.
    If this fails, the ``default`` keyword will be used.

    Parameters
    ----------
    url : `str`
        The file URL download.
    directory : `str`
        The directory path to download the file in to.
    default : `str`, optional
        The name to use if the first two methods fail. Defaults to "file".
    overwrite : `bool`, optional
        If `True` will overwrite a file of the same name. Defaults to `False`.

    Returns
    -------
    `str`:
        The file path for the downloaded file.
    """
    opn = urlopen(url)
    try:
        path = download_fileobj(opn, directory, url, default, overwrite)
    finally:
        opn.close()
    return path

# These two functions were in the stdlib cgi module which was deprecated in Python 3.11
# They are copied here under the terms of the PSF Licence 2.0


def _parseparam(s):
    while s[:1] == ';':
        s = s[1:]
        end = s.find(';')
        while end > 0 and (s.count('"', 0, end) - s.count('\\"', 0, end)) % 2:
            end = s.find(';', end + 1)
        if end < 0:
            end = len(s)
        f = s[:end]
        yield f.strip()
        s = s[end:]


def parse_header(line):
    """Parse a Content-type like header.

    Return the main content-type and a dictionary of options.

    """
    parts = _parseparam(';' + line)
    key = parts.__next__()
    pdict = {}
    for p in parts:
        i = p.find('=')
        if i >= 0:
            name = p[:i].strip().lower()
            value = p[i+1:].strip()
            if len(value) >= 2 and value[0] == value[-1] == '"':
                value = value[1:-1]
                value = value.replace('\\\\', '\\').replace('\\"', '"')
            pdict[name] = value
    return key, pdict
