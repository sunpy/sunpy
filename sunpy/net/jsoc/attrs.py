from __future__ import absolute_import
from sunpy.net.attr import (Attr, AttrWalker, AttrAnd, AttrOr )
from sunpy.net.vso.attrs import Time, _VSOSimpleAttr

__all__ = ['series', 'protocol','notify','compression','walker']
class series(_VSOSimpleAttr):
    pass

class protocol(_VSOSimpleAttr):
    pass

class notify(_VSOSimpleAttr):
    pass 

class compression(_VSOSimpleAttr):
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

    imap[query.__class__.__name__] = query.value

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

