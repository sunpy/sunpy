# -*- coding: utf-8 -*-
# Author: Florian Mayer <florian.mayer@bitsrc.org>

import os
import sys
import tempfile
import threading

from functools import partial
from collections import defaultdict

from suds import client
import sunpy
from sunpy.net import download

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

class Results(object):
    def __init__(self, callback, n=0, done=None):
        self.callback = callback
        self.n = n
        self.map_ = {}
        self.done = done
        self.evt = threading.Event()
    
    def submit(self, keys, value):
        for key in keys:
            self.map_[key] = value
        self.poke()
    
    def poke(self):
        self.n -= 1
        if not self.n:
            if self.done is not None:
                self.map_ = self.done(self.map_)
            self.callback(self.map_)
            self.evt.set()
    
    def require(self, keys):
        self.n += 1
        return partial(self.submit, keys)
    
    def wait(self):
        self.evt.wait()
        return self.map_


def mk_filename(pattern, response, sock, url):
    # FIXME: not name
    name = sock.headers.get(
        'Content-Disposition', url.rstrip('/').rsplit('/', 1)[-1]
    )
    name = pattern.format(file=name, **dict(response))
    dir_ = os.path.dirname(name)
    if not os.path.exists(dir_):
        os.makedirs(os.path.dirname(name))
    return name


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
    
    def query_legacy(self, time_start, time_end, **kwargs):
        ALIASES = {'wave_min': 'wave_wavemin', 'wave_max': 'wave_wavemax',
                   'wave_type': 'wave_wavetype', 'wave_unit': 'wave_waveunit'}
        kwargs.update({'time_start': time_start, 'time_end': time_end})
        
        queryreq = self.api.factory.create('QueryRequest')
        for k, v in kwargs.iteritems():
            k = ALIASES.get(k, k)
            if k.startswith('time'):
                v = v.strftime(TIMEFORMAT)
            attr = k.split('_')
            lst = attr[-1]
            rest = attr[:-1]
            
            item = queryreq.block
            for elem in rest:
                item = item[elem]
            item[lst] = v
        return self.api.service.Query(queryreq)
    
    def get(self, query_response, path=None, downloader=None, methods=['URL-FILE']):
        if downloader is None:
            downloader = download.Downloader(1)
            threading.Thread(target=downloader.reactor.run).start()
            res = Results(
                lambda _: downloader.reactor.stop(), 1,
                lambda mp: self.link(qr, mp)
            )
        else:
            res = Results(
                lambda _: None, 1, lambda mp: self.link(qr, mp)
            )
        if path is None:
            path = os.path.join(tempfile.mkdtemp(), '{file}')
        self.download_all(
            api.api.service.GetData(
                api.make_getdatarequest(query_response, methods)
                ),
            methods, downloader, path,
            API.by_fileid(query_response), res
        )
        res.poke()
        return res
    
    @staticmethod
    def link(query_response, map_):
        ret = []
        for prov_item in query_response.provideritem:
            for record_item in prov_item.record.recorditem:
                ret.append(Blah(record_item, map_[record_item.fileid]['path']))
        return ret
    
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
    
    def download_all(self, response, methods, dw, path, qr, res, info=None):
        for dresponse in response.getdataresponseitem:
            code = (
                dresponse.status[5:8] if hasattr(dresponse, 'status') else '200'
            )
            if code == '200':
                for dataitem in dresponse.getdataitem.dataitem:
                    location = self.download(
                        dresponse.method.methodtype[0],
                        dataitem.url,
                        dw,
                        res.require(map(str, dataitem.fileiditem.fileid)),
                        path,
                        qr[dataitem.fileiditem.fileid[0]]
                    )
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
                
                self.download_all(
                    self.api.service.GetData(request), methods, dw, path,
                    qr, res, info
                )
    
    def download(self, method, url, dw, callback, *args):
        if method.startswith('URL'):
            dw.reactor.call_sync(
                partial(dw.download, url, partial(mk_filename, *args),
                        callback)
            )
    
    @staticmethod
    def by_provider(response):
        map_ = defaultdict(list)
        for prov_item in response.provideritem:
            map_[prov_item.provider].extend(prov_item.record.recorditem)
        return map_
    
    @staticmethod
    def by_fileid(response):
        map_ = {}
        for prov_item in response.provideritem:
            for record_item in prov_item.record.recorditem:
                map_[record_item.fileid] = record_item
        return map_
    
    def multiple_choices(self, choices, response):
        for elem in self.method_order:
            if elem in choices:
                return [elem]
        raise ValueError
    
    def missing_information(self, info, field):
        print info, field
        return {'email': 'segfaulthunter@gmail.com'}


# TODO: class InteractiveAPI(API)

class Blah(object):
    def __init__(self, recorditem, path):
        self.recorditem = recorditem
        self.path = path
    
    def to_map(self):
        return sunpy.Map(self.path)


if __name__ == '__main__':
    from datetime import datetime
    api = API()
    
    qr = api.query_legacy(
        datetime(2010, 1, 1), datetime(2010, 1, 1, 1),
        instrument='eit'
    )
    
    qr = api.query(
        Time(datetime(2010, 1, 1), datetime(2010, 1, 1, 1)) & Instrument('eit')
    )
    res = api.get(qr, methods=['URL']).wait()
    res[0].to_map().plot()

