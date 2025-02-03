from sunpy.net._attrs import Detector, Time, Wavelength
from sunpy.net.attr import AttrAnd, AttrOr, AttrWalker, DataAttr, SimpleAttr


__all__ = ["Dataset", "Limit"]

class Dataset(SimpleAttr):
    """
    Name of the dataset required.
    """

class Limit(SimpleAttr):
    """
    Sets the max number of results, defaults to 20.
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


@walker.add_applier(Wavelength)
def _apply(wlk,query,imap):
    imap["wavelnth"] = int(query.min.value)


@walker.add_applier(Dataset)
def _apply(wlk,query,imap):
    imap["datasets"] = f"metadata_{query.value}"


@walker.add_applier(Time)
def _apply(wlk,query,imap):
    imap["date_end__gte"] = query.start.value
    imap['data_beg__lte'] = query.end.value


@walker.add_applier(Detector)
def _apply(wlk,query,imap):
    imap["detector__iexact"] = query.value


@walker.add_applier(Limit)
def _apply(wlk,query,imap):
    imap["limit"] = int(query.value)
