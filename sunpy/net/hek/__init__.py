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
from sunpy.util.xml import xml_to_dict

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
    attrs.walker.apply(attrs.SpatialRegion(), {}, default)
    
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
    @property
    def vso_time(self):
        return v_attrs.Time(
            datetime.strptime(self['event_starttime'], "%Y-%m-%dT%H:%M:%S"),
            datetime.strptime(self['event_endtime'], "%Y-%m-%dT%H:%M:%S")
        )
    
    @property
    def vso_instrument(self):
        if self['obs_instrument'] == 'HEK':
            raise ValueError("No instrument contained.")
        return v_attrs.Instrument(self['obs_instrument'])
    
    @property
    def vso_all(self):
        return attr.and_(self.vso_time, self.vso_instrument)
    
    def get_voevent(self, as_dict=True):
        """Retrieves the VOEvent object associated with a given event and
        returns it as either a Python dictionary or an XML string."""
        
        base_url = "http://www.lmsal.com/hek/her?"
        
        # Build URL
        params = {                                                      
            "cmd": "export-voevent",
            "cosec": 1,
            "ivorn": self['kb_archivid']
        }
        url = base_url + urlencode(params)
        
        # Query and read response
        response = urlopen(url).read()

        # Return a string or dict
        if as_dict:
            return xml_to_dict(response)
        else:
            return response


if __name__ == '__main__':
    import pprint
    from sunpy.net.hek import attrs as a

    c = HEKClient()
    b = c.query(
        a.Time((2010, 1, 1), (2010, 1, 2)) | a.Time((2010, 1, 3), (2010, 1, 4)),
        a.AR, a.FL
    )
    print b[0].vso_all
