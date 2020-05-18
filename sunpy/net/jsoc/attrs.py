from sunpy.net._attrs import Time, Wavelength
from sunpy.net.attr import AttrAnd, AttrOr, AttrWalker, DataAttr, SimpleAttr

__all__ = ['Series', 'Protocol', 'Notify', 'Segment', 'Keys', 'PrimeKey']


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
        return f"{object.__repr__(self)}" + "\n" + f"{self.label, self.value}"

    def collides(self, other):
        return False


class Segment(SimpleAttr):
    """
    Segments choose which files to download when there are more than
    one present for each record e.g. 'image'.
    """

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
