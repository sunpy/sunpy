from sunpy.net._attrs import Detector, Time, Wavelength
from sunpy.net.attr import AttrAnd, AttrWalker, DataAttr, SimpleAttr

__all__ = ["Dataset", "Limit" , "Target" , "Tags"]

class Dataset(SimpleAttr):
    """
    Name of the dataset required.
    """

class Limit(SimpleAttr):
    """
    Sets the max number of results, defaults to 20.
    """

class Target(SimpleAttr):
    """
    Indicates the observation region. Flare site used when flare flag is set.
    Source of information observation planning database, or telemetry if flare flag is set.
    some of the targets include
    Active Region (AH)
    Coronal Hole (CH)
    Flare (FS)
    Quiet Region (QR)
    """

class Tags(SimpleAttr):
    """
    a simple information associated to a solar observation,
    for example "moon transit"  etc.
    To each metadata instance can be associated 0 or more tags.
    """

walker = AttrWalker()

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

@walker.add_applier(Target)
def _apply(wlk,query,imap):
    imap["target__in"] = query.value

@walker.add_applier(Tags)
def _apply(wlk,query,imap):
    imap["tag"] = query.value
