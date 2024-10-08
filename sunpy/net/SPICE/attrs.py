import sunpy.net.attrs as a
from sunpy.net.attr import AttrAnd, AttrOr, AttrWalker, DataAttr, SimpleAttr


#PSP SPECIFIC
class Analysis_fk(SimpleAttr):
    """
    to get analysis frame kernel
    """
    def __init__(self,value):
        self.value = value
        if not isinstance(value,bool):
            raise ValueError("given value must be boolean")

class Numupdates(SimpleAttr):
    """
    Number of times file has been updated since launch
    """

#Solo SPECIFIC
class Voem(SimpleAttr):
    """
    Voem : reference to the source OEM file version
    """

class Readme(SimpleAttr):
    def __init__(self,value):

        if not isinstance(value,bool):
            raise ValueError(f"value must be boolean not {type(value)}")
        self.value = value

class Sensor(SimpleAttr):
    """
    sensor for kernels
    """

#general
class Kernel_type(SimpleAttr):
    """
    kernel type
    """
    def __init__(self, value):
        self.value = value
        if value is None:
            raise ValueError ("kernel type is required")

class Instrument(SimpleAttr):
    """
    Instrument for kernels
    """

class Link(SimpleAttr):
    """
    name of link for spice kernels
    """

class Version(SimpleAttr):
    """
    version number for kernels
    """


class Index(SimpleAttr):
    """
    index of link to be downloaded
    """
    def __init__(self,*value):
        self.value = value

class Obsevotory(SimpleAttr):
    def __init__(self, *value):
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

@walker.add_applier(Obsevotory)
def _(wlk,query,imap):
    imap[query.__class__.__name__.lower()] = query.value

@walker.add_applier(Kernel_type)
def _(wlk,query,imap):
    imap[query.__class__.__name__.lower()] = query.value

@walker.add_applier(a.Time)
def _(wlk,query,imap):
    imap["start"] = query.start.tdb.strftime('%Y%m%d')
    imap["end"] = query.end.tdb.strftime('%Y%m%d') if query.end is not None else None

@walker.add_applier(a.Instrument)
def _(wlk,query,imap):
    imap[query.__class__.__name__.lower()] = query.value

@walker.add_applier(AttrAnd)
def _apply(wlk, query, imap):
    for iattr in query.attrs:
        wlk.apply(iattr, imap)
