# -*- coding: utf-8 -*-
# Author: Florian Mayer <florian.mayer@bitsrc.org>

from __future__ import division
from __future__ import print_function
from __future__ import absolute_import
from __future__ import unicode_literals

import os
import re
import sys
import shutil

# For Content-Disposition parsing
from urllib2 import urlopen
from urlparse import urlparse
from email.parser import FeedParser
from unicodedata import normalize
from itertools import ifilter

from sunpy.util import replacement_filename

__all__ = ['slugify','get_content_disposition', 'get_filename',
           'get_system_filename', 'get_system_filename_slugify',
           'download_file', 'download_fileobj']

# Characters not allowed in slugified version.
_punct_re = re.compile(r'[:\t !"#$%&\'()*\-/<=>?@\[\\\]^_`{|},.]+')

def slugify(text, delim=u'_', encoding="ascii"):
    """ Slugify given unicode text. """
    text = normalize('NFKD', text)
    return unicode(delim).join(ifilter(None, (
        word.encode(encoding, 'ignore')
        for word in _punct_re.split(text.lower())        
        )))


def get_content_disposition(content_disposition):
    """ Get content disposition filename from given header. Do not include
    "Content-Disposition:". Returns a unicode string! """
    parser = FeedParser()
    parser.feed(b'Content-Disposition: ' + content_disposition)
    name = parser.close().get_filename()
    if not isinstance(name, unicode):
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
    return unicode(name)


def get_system_filename(sock, url, default=u"file"):
    """ Get filename from given urllib2.urlopen object and URL.
    First, attempts to extract Content-Disposition, second, extract
    from URL, eventually fall back to default. Returns bytestring
    in file system encoding. """
    name = get_filename(sock, url)
    if not name:
        name = default.decode("ascii", "ignore")
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


def download_fileobj(opn, directory, url='', default=u"file", overwrite=False):
    """ Download file from url into directory. Try to get filename from
    Content-Disposition header, otherwise get from path of url if given.
    Fall back to default if both fail. Only overwrite existing files when
    overwrite is True. """
    filename = get_system_filename(opn, url, default)
    path = os.path.join(directory, filename)
    if not overwrite and os.path.exists(path):
        path = replacement_filename(path)
    with open(path, 'wb') as fd:
        shutil.copyfileobj(opn, fd)
    return path
