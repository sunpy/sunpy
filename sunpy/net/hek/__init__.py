# -*- coding: utf-8 -*-
# Author: Florian Mayer <florian.mayer@bitsrc.org>

# pylint: disable=C0103,R0903

from __future__ import absolute_import

import json

from itertools import chain
from urllib2 import urlopen
from urllib import urlencode
from datetime import datetime
from functools import partial

from sunpy.net import attr
from sunpy.net.hek import attrs
from sunpy.net.vso import attrs as v_attrs
from sunpy.util.util import unique

DEFAULT_URL = 'http://www.lmsal.com/hek/her'

def _freeze(obj):
    """ Create hashable representation of result dict. """
    if isinstance(obj, dict):
        return tuple((k, _freeze(v)) for k, v in obj.iteritems())
    if isinstance(obj, list):
        return tuple(_freeze(elem) for elem in obj)
    return obj


class HEKClient(object):
    # FIXME: Types!
    
    default = {
        'cosec': '2',
        'cmd': 'search',
        'type': 'column',
        'event_type': '**',
    }
    # Default to full disk.
    attrs.walker.apply(attrs.SpartialRegion(), {}, default)
    
    def __init__(self, url=DEFAULT_URL):
        self.url = url
    
    def _download(self, data):
        page = 1        
        results = []
        
        while True:
            data['page'] = page
            result = json.load(urlopen(self.url, urlencode(data)))
            results.extend(result['result'])
            
            if not result['overmax']:
                return map(Response, results)
            page += 1
    
    def query(self, *query):
        query = attr.and_(*query)
        
        data = attrs.walker.create(query, {})
        ndata = []
        for elem in data:
            new = self.default.copy()
            new.update(elem)
            ndata.append(new)
        
        if len(ndata) == 1:
            return self._download(ndata[0])
        else:
            return self.merge(self._download(data) for data in ndata)
    
    def merge(self, responses):
        return list(unique(chain.from_iterable(responses), _freeze))


class Response(dict):
    def to_vso(self):
        return v_attrs.Time(
            datetime.strptime(self['event_starttime'], "%Y-%m-%dT%H:%M:%S"),
            datetime.strptime(self['event_endtime'], "%Y-%m-%dT%H:%M:%S")
        )


if __name__ == '__main__':
    import pprint
    from sunpy.net.hek import attrs as a

    c = HEKClient()
    print len(c.query(
        a.Time((2010, 1, 1), (2010, 1, 2)) | a.Time((2010, 1, 3), (2010, 1, 4)),
        a.AR, a.FL
    ))
