# -*- coding: utf-8 -*-
# Author: Florian Mayer <florian.mayer@bitsrc.org>

import sys
import os
from collections import defaultdict

from suds import client

DEFAULT_URL = 'http://docs.virtualsolar.org/WSDL/VSOi_rpc_literal.wsdl'
TIMEFORMAT = '%Y%m%d%H%M%S'


from suds.transport import Reply
from suds.transport.https import HttpAuthenticated

class LoopbackTransport(HttpAuthenticated):
    def send(self, request):
        return Reply(200, {}, request.message)


loopback = client.Client(DEFAULT_URL, transport=LoopbackTransport(), retxml=True)


class _Attr(object):
    def __and__(self, other):
        return _AttrAnd([self, other])


class _AttrAnd(_Attr):
    def __init__(self, attrs):
        self.attrs = attrs
    
    def __and__(self, other):
        return _AttrAnd(self.attrs + [other])
    
    def apply(self, queryblock):
        for attr in self.attrs:
            attr.apply(queryblock)


class _ComplexAttr(_Attr):
    def __init__(self, attr, attrs):
        self.attr = attr
        self.attrs = attrs
    
    def apply(self, queryblock):
        for name in self.attr:
            queryblock = getattr(queryblock, name)
        
        for k, v in self.attrs.iteritems():
            queryblock[k] = v


class _SimpleAttr(_Attr):
    def __init__(self, field, value):
        self.field = field
        self.value = value
    
    def apply(self, queryblock):
        queryblock[self.field] = self.value

# ----------------------------------------

class Wave(_ComplexAttr):
    def __init__(self, wavemin, wavemax, waveunit):
        _ComplexAttr.__init__(self, ['wave'], {
            'wavemin': wavemin,
            'wavemax': wavemax,
            'waveunit': waveunit,
        })


class Time(_ComplexAttr):
    def __init__(self, start, end, near=None):
        _ComplexAttr.__init__(self, ['time'], {
            'start': start.strftime(TIMEFORMAT),
            'end': end.strftime(TIMEFORMAT),
            'near': near.strftime(TIMEFORMAT) if near is not None else '',
        })


class Extent(_ComplexAttr):
    def __init__(self, x, y, width, length, type_):
        _ComplexAttr.__init__(self, ['extent'], {
            'x': x,
            'y': y,
            'width': width,
            'length': length,
            'type': type_,
        })


class Field(_ComplexAttr):
    def __init__(self, fielditem):
        _ComplexAttr.__init__(self, ['field'], {
            'fielditem': fielditem
        })


def  _mk_simpleattr(field):
    """ Create a _SimpleField class for field. """
    class _foo(_SimpleAttr):
        def __init__(self, arg):
            _SimpleAttr.__init__(self, field, arg)
    _foo.__name__ = field.capitalize()
    return _foo


for elem in ['provider', 'source', 'instrument', 'physobs', 'pixels',
             'level', 'resolution', 'detector', 'filter', 'sample',
             'quicklook', 'pscale']:
    setattr(sys.modules[__name__], elem.capitalize(), _mk_simpleattr(elem))

# ----------------------------------------


class API(object):
    method_order = [
        'URL-TAR_GZ', 'URL-ZIP', 'URL-TAR', 'URL-FILE', 'URL-packaged'
    ]
    def __init__(self, url=DEFAULT_URL, *args):
        self.api = client.Client(url, *args)
        self.api.set_options(port='nsoVSOi')
    
    def query(self, query):
        queryreq = self.api.factory.create('QueryRequest')
        query.apply(queryreq.block)
        return self.api.service.Query(queryreq)
    
    def make_getdatarequest(self, response, methods=None):
        if methods is None:
            methods = self.method_order + ['URL']
        
        return self.create_getdatarequest(
            dict((k, [x.fileid for x in v])
                 for k, v in self.by_provider(response).iteritems()),
            methods
        )
    
    def create_getdatarequest(self, map_, methods, info=None):
        if info is None:
            info = {}
        
        request = self.api.factory.create('VSOGetDataRequest')
        r = request.request
        r.method.methodtype.extend(methods)
        
        for k, v in info.iteritems():
            r.info[k] = v
        
        for k, v in map_.iteritems():
            datarequest = self.api.factory.create('DataRequestItem')
            datarequest.provider = k
            datarequest.fileiditem.fileid = v
            r.datacontainer.datarequestitem.append(datarequest)
        
        return request
    
    def download_all(self, response, methods, info=None):
        map_ = {}
        for dresponse in response.getdataresponseitem:
            code = (
                dresponse.status[5:8] if hasattr(dresponse, 'status') else '200'
            )
            if code == '200':
                for dataitem in dresponse.getdataitem.dataitem:
                    path = self.download(
                        dresponse.method.methodtype[0],
                        dataitem.url
                    )
                    for fileid in dataitem.fileiditem.fileid:
                        map_[fileid] = path
            elif code == '300' or code == '412':
                files = []
                for dataitem in dresponse.getdataitem.dataitem:
                    files.extend(dataitem.fileiditem.fileid)
                if code == '300':
                    try:
                        methods = self.multiple_choices(
                            dresponse.method.methodtype, dresponse
                        )
                    except ValueError:
                        # TODO: Log.
                        continue
                elif code == '412':
                    info = self.missing_information(
                        info, dresponse.info
                    )
                request = self.create_getdatarequest(
                    {dresponse.provider: files}, methods, info
                )
                
                map_.update(
                    self.download_all(
                        self.api.service.GetData(request), methods, info
                    )
                )
        return map_
    
    def download(self, method, url):
        print method, url
    
    @staticmethod
    def by_provider(response):
        map_ = defaultdict(list)
        for prov_item in a.provideritem:
            map_[prov_item.provider].extend(prov_item.record.recorditem)
        return map_
    
    def multiple_choices(self, choices, response):
        for elem in self.method_order:
            if elem in choices:
                return [elem]
        raise ValueError
    
    def missing_information(self, info, field):
        print info, field
        return {'EMAIL': 'segfaulthunter@gmail.com'}


# TODO: class InteractiveAPI(API)

if __name__ == '__main__':
    from datetime import datetime
    api = API()
    
    a = api.query(
        Time(datetime(2010, 1, 1), datetime(2010, 1, 1, 1)) & Instrument('eit')
    )
    
    print api.download_all(
        api.api.service.GetData(api.make_getdatarequest(a, ['URL'])),
        ['STAGING-ZIP']
    )
