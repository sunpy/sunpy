from __future__ import absolute_import

import numpy as np
import astropy.units as u

from sunpy.net.attr import (Attr, AttrWalker, AttrAnd, AttrOr)
from sunpy.net.vso.attrs import Time, _VSOSimpleAttr

__all__ = ['Series', 'Protocol', 'Notify', 'Compression', 'Wavelength', 'Time',
           'Segment', 'walker']


class Time(Time):
    """
    Time range to download
    """


class Series(_VSOSimpleAttr):
    """
    The JSOC Series to Download.

    See `this<http://jsoc.stanford.edu/JsocSeries_DataProducts_map.html>_`
    for a list of series'.
    """
    pass


class Segment(_VSOSimpleAttr):
    """
    Segments choose which files to download when there are more than
    one present for each record e.g. 'image'
    """
    pass


class Protocol(_VSOSimpleAttr):
    """
    The type of download to request one of
    ("FITS", "JPEG", "MPG", "MP4", or "as-is").
    Only FITS is supported, the others will require extra keywords.
    """
    pass


class Notify(_VSOSimpleAttr):
    """
    An email address to get a notification to when JSOC has staged your request
    """
    pass


class Compression(_VSOSimpleAttr):
    """
    Compression format for requested files.

    'rice' or None, download FITS files with RICE compression.
    """
    pass


class Wavelength(_VSOSimpleAttr):
    """
    Wavelength or list of wavelengths to download. Must be specified in correct
    units for the series.
    """
    def __init__(self, value):
        if not isinstance(value, u.Quantity):
            raise TypeError("Wave inputs must be astropy Quantities")
        Attr.__init__(self)

        self.value = int(np.ceil(value.to(u.AA).value))

    def __or__(self, other):
        if isinstance(other, self.__class__):
            return self.__class__([self.value, other.value])
        if self == other:
            return self
        return AttrOr([self, other])
    __ror__ = __or__


walker = AttrWalker()


@walker.add_creator(AttrAnd, _VSOSimpleAttr, Time)
def _create(wlk, query):

    map_ = {}
    wlk.apply(query, map_)
    return [map_]


@walker.add_applier(AttrAnd)
def _apply(wlk, query, imap):

    for iattr in query.attrs:
        wlk.apply(iattr, imap)


@walker.add_applier(_VSOSimpleAttr)
def _apply(wlk, query, imap):

    imap[query.__class__.__name__.lower()] = query.value


@walker.add_applier(Time)
def _apply(wlk, query, imap):

    imap['start_time'] = query.start
    imap['end_time'] = query.end


@walker.add_creator(AttrOr)
def _create(wlk, query):

    qblocks = []
    for iattr in query.attrs:
        qblocks.extend(wlk.create(iattr))

    return qblocks
