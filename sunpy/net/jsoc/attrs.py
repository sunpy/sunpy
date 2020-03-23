# -*- coding: utf-8 -*-
from astropy.time import Time as astropyTime

from sunpy.net._attrs import Wavelength, Time
from sunpy.net.attr import AttrWalker, AttrAnd, AttrOr, DataAttr, SimpleAttr


__all__ = ['Series', 'Protocol', 'Notify', 'Segment', 'Keys', 'PrimeKey', 'Quality']


class Series(SimpleAttr):
    """
    The JSOC Series to Download.

    This is the list of `Series <http://jsoc.stanford.edu/JsocSeries_DataProducts_map.html>`__.
    """


class Keys(SimpleAttr):
    """
    Keys choose which keywords to fetch while making a query request.
    """


class PrimeKey(DataAttr):
    """
    Prime Keys

    Parameters
    ----------
    label : str
    value : str
    """
    def __init__(self, label, value):
        super().__init__()
        self.label = label
        self.value = value

    def __repr__(self):
        return "<{cname!s}({lab!r},{val!r})>".format(
            cname=self.__class__.__name__, lab=self.label, val=self.value)

    def collides(self, other):
        return False


class Segment(SimpleAttr):
    """
    Segments choose which files to download when there are more than
    one present for each record e.g. 'image'.
    """
    def __init__(self, value):
        super().__init__(value)
        self.value = value

    def __repr__(self):
        return "<{cname!s}({val!r})>".format(
            cname=self.__class__.__name__, val=self.value)

    def collides(self, other):
        return False


class Protocol(SimpleAttr):
    """
    The type of download to request one of
    ("FITS", "JPEG", "MPG", "MP4", or "as-is").
    Only FITS is supported, the others will require extra keywords.
    """


class Notify(SimpleAttr):
    """
    An email address to get a notification to when JSOC has staged your request.
    """
    def __init__(self, value):
        super().__init__(value)
        if value.find('@') == -1:
            raise ValueError("Notify attribute must contain an '@' symbol "
                             "to be a valid email address")
        self.value = value


class Quality(SimpleAttr):
    """
    A comparison test and number to run against the quality value (e.g. `>= 0`).

    Quality information is encoded into 32 bit flags e.g. `0b01010101010101010101010101010101`
    which can be represented in hexadecimal as `0x55555555`,  or as in integer `1431655765`.
    The exact meaning of each bit defers from series to series see
    `MDI <http://jsoc.stanford.edu/MDI/MDI_Global.html>`_,
    `HMI <http://jsoc.stanford.edu/doc/data/hmi/Quality_Bits/quallev1.txt>`_, and
    `AIA <http://jsoc.stanford.edu/~jsoc/keywords/AIA/AIA02840_K_AIA-SDO_FITS_Keyword_Document.pdf>`_.

    .. warning::
        The quality field is returned as an **unsigned 32 bit** integer representation of the bit
        flags but the attribute expects a **signed 32 bit**  integer representation of the bit
        flags.

    For example the bit flags `0b10000000000000000000000000000000` or `0x80000000` in hex have an
    unsigned 32 bit integer value of `2147483648` or a signed 32 bit integer value of `-2147483648`.

    For HMI the missing data flag is `0x80000000` or `0b10000000000000000000000000000000` this can
    be converted to a signed int using::

        import numpy as np
        np.array(0x80000000).astype(np.int32)
        array(-2147483648, dtype=int32)

    This value could then be used with the quality attribute `Quality('>-2147483648)` however this
    may not remove all missing data as other flag or bits may be set. For singed 32 bit these
    numbers would be smaller (`0x80000000 < 0x88000000`, `-2147483648 < -2013265920`) and so the
    comparision test would be false. In this case it is important to note all signed 32 bit numbers
    with the top bit set (`0x80000000`) are negative so it is sufficient to query for quality values
    the are greater than or equal to zero `Quality('>=0')`

    Another example with HMI the QUAL_ECLIPSE flag (`0x00000200`) as a signed 32 bit number is `512`
    setting any additional flags can only increase this value so to filter out any results with the
    QUAL_ECLIPSE flag `Quality( < 512)` would be used.
    """


walker = AttrWalker()


@walker.add_creator(AttrOr)
def _create1(wlk, query):

    qblocks = []
    for iattr in query.attrs:
        qblocks.extend(wlk.create(iattr))

    return qblocks


@walker.add_creator(AttrAnd, DataAttr)
def _create(wlk, query):

    map_ = {}
    wlk.apply(query, map_)
    return [map_]


@walker.add_applier(AttrAnd)
def _apply(wlk, query, imap):

    for iattr in query.attrs:
        wlk.apply(iattr, imap)


@walker.add_applier(SimpleAttr)
def _apply1(wlk, query, imap):

    imap[query.__class__.__name__.lower()] = query.value


@walker.add_applier(PrimeKey)
def _apply1(wlk, query, imap):

    key = 'primekey'
    if key in imap:
        imap[key][query.label] = query.value
    else:
        imap[key] = {query.label: query.value}


@walker.add_applier(Segment)
def _apply1(wlk, query, imap):

    key = 'segment'
    if key in imap:
        imap[key].append(query.value)
    else:
        imap[key] = [query.value]


@walker.add_applier(Time)
def _apply2(wlk, query, imap):
    imap['start_time'] = query.start
    imap['end_time'] = query.end


@walker.add_applier(Wavelength)
def _apply_wave(wlk, query, imap):
    if query.min != query.max:
        raise ValueError(
            "For JSOC queries Wavelength.min must equal Wavelength.max")

    imap[query.__class__.__name__.lower()] = query.min
