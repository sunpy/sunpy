from __future__ import absolute_import

__all__ = ['Series', 'Protocol', 'Notify', 'Compression', 'Wavelength', 'Time',
           'Segment', 'Sample']

import astropy.units as u

from sunpy.net.attr import Attr, AttrWalker, AttrAnd, AttrOr
from sunpy.net.vso.attrs import _VSOSimpleAttr
from sunpy.net.vso.attrs import Time as vTime, Sample as vSample

###############################################################################
# This is a horrific hack to make automodapi pick up these as jsoc attrs.


class Time(vTime):
    __doc__ = vTime.__doc__
    pass


class Sample(vSample):
    __doc__ = vSample.__doc__
    pass

###############################################################################


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
    def __init__(self, value):
        super(Notify, self).__init__(value)
        if value.find('@') == -1:
            raise ValueError("Notify attribute must contain an '@' symbol "
                             "to be a valid email address")
        self.value = value


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
        if not (isinstance(value, u.Quantity) or isinstance(value, list)):
            raise TypeError("Wave inputs must be astropy Quantities")
        Attr.__init__(self)

        self.value = value

    def __or__(self, other):
        if self == other:
            return self
        if isinstance(other, self.__class__):
            return self.__class__([self.value, other.value])
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
def _apply1(wlk, query, imap):

    imap[query.__class__.__name__.lower()] = query.value


@walker.add_applier(Time)
def _apply2(wlk, query, imap):

    imap['start_time'] = query.start
    imap['end_time'] = query.end


@walker.add_creator(AttrOr)
def _create1(wlk, query):

    qblocks = []
    for iattr in query.attrs:
        qblocks.extend(wlk.create(iattr))

    return qblocks
