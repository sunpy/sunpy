"""
This module provides general net utility functions.
"""
import os
import re
import sys
import shutil
from unicodedata import normalize
from email.parser import FeedParser
from urllib.parse import urljoin, urlparse
from urllib.request import urlopen

from sunpy.util import replacement_filename

__all__ = ['slugify', 'get_content_disposition', 'get_filename', 'get_system_filename',
           'download_file', 'download_fileobj', 'check_download_file']

# Characters not allowed in slugified version.
_punct_re = re.compile(r'[:\t !"#$%&\'()*\-/<=>?@\[\\\]^_`{|},.]+')


def slugify(text, delim='_'):
    """
    Slugify given unicode text.

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
    text = normalize('NFKD', text)

    period = '.'

    name_and_extension = text.rsplit(period, 1)
    name = name_and_extension[0]

    name = str(delim).join(
        filter(None, (word for word in _punct_re.split(name.lower()))))

    if len(name_and_extension) == 2:
        extension = name_and_extension[1]
        return str(period).join([name, extension])
    else:
        return name


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
    `str`:
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
    overwrite: `bool`, optional
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
    overwrite: `bool`, optional
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


def check_download_file(filename, remotepath, download_dir, remotename=None, replace=False):
    """
    Downloads a file from a remotepath to a localpath if it isn't there.

    This function checks whether a file with name ``filename`` exists in the user's local machine.
    If it doesn't, it downloads the file from ``remotepath``.

    Parameters
    ----------
    filename : `str`
        Name of file.
    remotepath : `str`
        URL of the remote location from which filename can be downloaded.
    download_dir : `str`
        The files directory.
    remotename : `str`, optional
        Filename under which the file is stored remotely.
        Default is same as filename.
    replace : `bool`, optional
        If `True`, file will be downloaded whether or not file already exists locally.

    Examples
    --------
    >>> from sunpy.util.net import check_download_file
    >>> remotepath = "https://www.download_repository.com/downloads/"
    >>> check_download_file("filename.txt", remotepath, download_dir='.')  # doctest: +SKIP
    """
    # Check if file already exists locally.  If not, try downloading it.
    if replace or not os.path.isfile(os.path.join(download_dir, filename)):
        # set local and remote file names be the same unless specified
        # by user.
        if not isinstance(remotename, str):
            remotename = filename

        download_file(urljoin(remotepath, remotename),
                      download_dir, default=filename, overwrite=replace)
