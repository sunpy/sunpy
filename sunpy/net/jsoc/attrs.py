from __future__ import absolute_import
from sunpy.net.attr import (Attr, AttrWalker, AttrAnd, AttrOr )
from sunpy.net.vso.attrs import Time, _VSOSimpleAttr

__all__ = ['Series', 'Protocol','Notify','Compression','walker']
class Series(_VSOSimpleAttr):
    pass

class Protocol(_VSOSimpleAttr):
    pass

class Notify(_VSOSimpleAttr):
    pass 

class Compression(_VSOSimpleAttr):
    pass

walker = AttrWalker()

@walker.add_creator(AttrAnd,_VSOSimpleAttr,Time)
def _create(wlk,query):

    map_ = {}
    wlk.apply(query,map_)
    return [map_]

@walker.add_applier(AttrAnd)
def _apply(wlk,query,imap):
    
    for iattr in query.attrs:
        wlk.apply(iattr,imap)

@walker.add_applier(_VSOSimpleAttr)
def _apply(wlk,query,imap):

    imap[query.__class__.__name__.lower()] = query.value

@walker.add_applier(Time)
def _apply(wlk,query,imap):

    imap['start_time'] = query.start
    imap['end_time'] = query.end

     
@walker.add_creator(AttrOr)
def _create(wlk,query):

   qblocks = []
   for iattr in query.attrs:
       qblocks.extend(wlk.create(iattr))

   return qblocks

