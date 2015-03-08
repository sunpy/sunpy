from __future__ import absolute_import

import warnings

import astropy.units as u
from astropy.tests.helper import assert_quantity_allclose

from sunpy.net.attr import AttrWalker, AttrAnd, AttrOr
from sunpy.net.vso.attrs import Time, _VSOSimpleAttr, Wavelength

__all__ = ['Series', 'Protocol', 'Notify', 'Compression', 'Wavelength', 'Time',
           'Segment', 'walker']

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
            raise ValueError("Notify attribute must contain an '@' symbol to be a valid email address")

        self.value = value
    pass


class Compression(_VSOSimpleAttr):
    """
    Compression format for requested files.

    'rice' or None, download FITS files with RICE compression.
    """
    pass

walker = AttrWalker()


@walker.add_creator(AttrAnd, _VSOSimpleAttr, Time, Wavelength)
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

@walker.add_applier(Wavelength)
def _apply(wlk, query, imap):
    if assert_quantity_allclose(query.min, query.max):
        imap['wavelength'] = query.min.to(u.AA, equivalencies=u.spectral()).value
    else:
        raise ValueError('Minimum and Maximum wavelength have to be identical for JSOC queries.')


@walker.add_creator(AttrOr)
def _create(wlk, query):

    qblocks = []
    for iattr in query.attrs:
        qblocks.extend(wlk.create(iattr))

    return qblocks
