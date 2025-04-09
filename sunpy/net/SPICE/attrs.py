import sunpy.net.attrs as a
from sunpy.net.attr import AttrAnd, AttrOr, AttrWalker, SimpleAttr


class Link(SimpleAttr):
    """
    name of link for spice kernels
    """

class Observatory(SimpleAttr):
    """
    Name of the mission
    """

class Leapseconds(SimpleAttr):
    """
    pass
    """
    def __init__(self, value = True):
        super().__init__(value)


walker = AttrWalker()

@walker.add_creator(AttrOr)
def _create1(wlk, query):
    qblocks = []
    for iattr in query.attrs:
        qblocks.extend(wlk.create(iattr))
    return qblocks


@walker.add_creator(AttrAnd, SimpleAttr)
def _create(wlk, query):
    map_ = {}
    wlk.apply(query, map_)
    return [map_]

@walker.add_applier(a.Provider)
def _(wlk,query,imap):
    imap[query.__class__.__name__.lower()] = query.value

@walker.add_converter(a.Time)
def _(wlk,query,imap):
    imap["start"] = query.start.tdb.strftime('%Y%m%d')
    imap["end"] = query.end.tdb.strftime('%Y%m%d') if query.end is not None else ""

@walker.add_applier(a.Instrument)
def _(wlk,query,imap):
    imap[query.__class__.__name__.lower()] = query.value

@walker.add_applier(Observatory)
def _(wlk,query,imap):
    imap[query.__class__.__name__.lower()] = query.value

@walker.add_applier(Link)
def _(wlk,query,imap):
    imap[query.__class__.__name__.lower()] = query.value

@walker.add_applier(AttrAnd)
def _apply(wlk, query, imap):
    for iattr in query.attrs:
        wlk.apply(iattr, imap)

@walker.add_applier(Leapseconds)
def _apply(wlk,query,imap):
    imap[query.__class__.__name__.lower()] = query.value