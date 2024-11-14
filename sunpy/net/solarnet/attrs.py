from sunpy.net._attrs import Detector, Time, Wavelength
from sunpy.net.attr import AttrAnd, AttrOr, AttrWalker, DataAttr, SimpleAttr


class Datasets(SimpleAttr):
    """
    Name of datasets
    """

class limit(SimpleAttr):
    """
    number of links needed (default value is 10)
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

@walker.add_applier(Datasets)
def _apply(wlk,query,imap):
    imap["datasets"] = f"metadata_{query.value}"

@walker.add_applier(Time)
def _apply(wlk,query,imap):
    if query.end is not None:
        imap["date_obs__gte"] = query.start.value
        imap['data_obs_lt'] = query.end.value
    else:
        imap["date_obs"] = query.start.value

@walker.add_applier(Detector)
def _apply(wlk,query,imap):
    imap["detector__iexact"] = query.value

@walker.add_applier(limit)
def _apply(wlk,query,imap):
    imap["limit"] = query.value
