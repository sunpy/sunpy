# -*- coding: utf-8 -*-
# Author: Florian Mayer <florian.mayer@bitsrc.org>

import os
import sys
import tempfile
import threading

from datetime import datetime
from functools import partial
from collections import defaultdict

from suds import client
from sunpy.net import download

DEFAULT_URL = 'http://docs.virtualsolar.org/WSDL/VSOi_rpc_literal.wsdl'
DEFAULT_PORT = 'nsoVSOi'
TIMEFORMAT = '%Y%m%d%H%M%S'


class _Str(str):
    pass


class _Attr(object):
    def __and__(self, other):
        if isinstance(other, _AttrOr):
            return _AttrOr([elem & self for elem in other.attrs])
        return _AttrAnd([self, other])
    
    def __or__(self, other):
        return _AttrOr([self, other])


class _DummyAttr(_Attr):
    def __and__(self, other):
        return other
    
    def __or__(self, other):
        return other
    
    def apply(self, queryblock):
        pass
    
    def create(self, api):
        return api.factory.create('QueryRequestBlock')


class _AttrAnd(_Attr):
    def __init__(self, attrs):
        self.attrs = attrs
    
    def __and__(self, other):
        if isinstance(other, _AttrAnd):
            return _AttrAnd(self.attrs + other.attrs)
        if isinstance(other, _AttrOr):
            return _AttrOr([elem & self for elem in other.attrs])
        return _AttrAnd(self.attrs + [other])
    
    __rand__ = __and__
    
    def apply(self, queryblock):
        for attr in self.attrs:
            attr.apply(queryblock)
    
    def create(self, api):
        # TODO: Prove that we can assume that only _SimpleAttr and
        # _ComplexAttr can possibly exist here.
        value = api.factory.create('QueryRequestBlock')
        self.apply(value)
        return [value]
    
    def __repr__(self):
        return "<_AttrAnd(%r)>" % self.attrs


class _AttrOr(_Attr):
    def __init__(self, attrs):
        self.attrs = attrs
    
    def __or__(self, other):
        if isinstance(other, _AttrOr):
            return _AttrOr(self.attrs + other.attrs)
        return _AttrOr(self.attrs + [other])
    
    __ror__ = __or__
    
    def __and__(self, other):
        return _AttrOr([elem & other for elem in self.attrs])
    
    __rand__ = __and__
    
    def __xor__(self, other):
        new = _DummyAttr()
        for elem in self.attrs:
            new |= elem ^ other
        return new
    
    def apply(self, queryblock):
        # TODO: Prove this is unreachable.
        raise NotImplementedError
    
    def create(self, api):
        blocks = []
        for attr in self.attrs:
            blocks.extend(attr.create(api))
        return blocks
    
    def __repr__(self):
        return "<_AttrOr(%r)>" % self.attrs


class _ComplexAttr(_Attr):
    def __init__(self, attr, attrs):
        self.attr = attr
        self.attrs = attrs
    
    def apply(self, queryblock):
        for name in self.attr:
            queryblock = getattr(queryblock, name)
        
        for k, v in self.attrs.iteritems():
            queryblock[k] = v
    
    def create(self, api):
        value = factory.create('QueryRequestBlock')
        self.apply(value)
        return [value]
    
    def __repr__(self):
        return "<_ComplexAttr(%r, %r)>" % (self.attr, self.attrs)


class _SimpleAttr(_Attr):
    def __init__(self, field, value):
        self.field = field
        self.value = value
    
    def apply(self, queryblock):
        queryblock[self.field] = self.value
    
    def create(self, api):
        value = factory.create('QueryRequestBlock')
        self.apply(value)
        return value
    
    def __repr__(self):
        return "<_SimpleAttr(%r, %r)>" % (self.field, self.value)

# ----------------------------------------

class Wave(_ComplexAttr):
    wavelength = [
        ('angstrom', 1e-10),
        ('nm', 1e-9),
        ('micron', 1e-6),
        ('mm', 1e-3),
        ('cm', 1e-2),
        ('m', 1e-6),
    ]
    energy = [
        ('ev', 1),
        ('kev', 1e3),
        ('mev', 1e6),
    ]
    frequency = [
        ('hz', 1),
        ('khz', 1e3),
        ('mhz', 1e6),
        ('ghz', 1e9),
    ]
    units = {}
    for k, v in wavelength:
        units[k] = ('wavelength', v)
    for k, v in energy:
        units[k] = ('energy', v)
    for k, v in frequency:
        units[k] = ('frequency', v)
    
    def __init__(self, wavemin, wavemax, waveunit):
        wavemin, wavemax = self.to_angstrom(wavemin, wavemax, waveunit)
        waveunit = 'Angstrom'
        
        self.min = wavemin
        self.max = wavemax
        self.unit = waveunit
        
        _ComplexAttr.__init__(self, ['wave'], {
            'wavemin': wavemin,
            'wavemax': wavemax,
            'waveunit': waveunit,
        })
    
    def __xor__(self, other):
        if self.unit != other.unit:
            return NotImplemented
        new = _DummyAttr()
        if self.min < other.min:
            new |= Wave(self.min, min(other.min, self.max), self.unit)
        if other.max < self.max:
            new |= Wave(other.max, self.max, self.unit)
        return new
    
    @classmethod
    def to_angstrom(cls, min_, max_, unit):
        # Speed of light in m/s.
        C = 299792458
        ANGSTROM = cls.units['angstrom'][1]
        
        try:
            type_, n = cls.units[unit.lower()]
        except KeyError:
            raise ValueError('Cannot convert %s to Angstrom' % unit)
        
        if type_ == 'wavelength':
            k = n / ANGSTROM
            return (min_ / k, self.max_ / k)
        if type_ == 'frequency':
            k = n / ANGSTROM
            return k * (C / max_), k * (C / min_)
        if type_ == 'energy':
            k = n / (ANGSTROM / 1e-2)
            return k * (1 / (8065.53 * max_)), k * (1 / (8065.53 * min_))
        else:
            raise ValueError('Unable to convert %s to Angstrom' % type_)


class Time(_ComplexAttr):
    def __init__(self, start, end, near=None):
        self.start = start
        self.end = end
        self.near = near
        
        _ComplexAttr.__init__(self, ['time'], {
            'start': start.strftime(TIMEFORMAT),
            'end': end.strftime(TIMEFORMAT),
            'near': near.strftime(TIMEFORMAT) if near is not None else '',
        })
    
    def __xor__(self, other):
        new = _DummyAttr()
        if self.start < other.start:            
            new |= Time(self.start, min(other.start, self.end))
        if other.end < self.end:
            new |= Time(other.end, self.end)
        return new
    
    @classmethod
    def dt(cls, start, end, near=None):
        if near is not None:
            near = datetime(*near)
        return cls(datetime(*start), datetime(*end), near)


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
    # FIXME: os.path.exists(name)
    name = sock.headers.get(
        'Content-Disposition', url.rstrip('/').rsplit('/', 1)[-1]
    )
    if not name:
        name = response.fileid.replace('/', '_')
    
    fname = pattern.format(file=name, **dict(response))
    dir_ = os.path.dirname(fname)
    if not os.path.exists(dir_):
        os.makedirs(os.path.dirname(fname))
    return fname


def _mk_queryreq(api, block):
    queryreq = api.factory.create('QueryRequest')
    queryreq.block = block
    return queryreq


class API(object):
    method_order = [
        'URL-TAR_GZ', 'URL-ZIP', 'URL-TAR', 'URL-FILE', 'URL-packaged'
    ]
    def __init__(self, api=None):
        if api is None:
            api = client.Client(DEFAULT_URL)
            api.set_options(port=DEFAULT_PORT)
        self.api = api
    
    def query(self, query):
        queryreq = self.api.factory.create('QueryRequest')
        query.create(self.api)
        return self.merge(
            self.api.service.Query(_mk_queryreq(self.api, block))
            for block in query.create(self.api)
        )
    
    def merge(self, queryresponses):      
        fileids = set()
        providers = {}
        
        for queryresponse in queryresponses:
            for provideritem in queryresponse.provideritem:
                provider = provideritem.provider
                if not provideritem.provider in providers:
                    providers[provider] = provideritem
                    fileids |= set(
                        record_item.fileid
                        for record_item in provideritem.record.recorditem
                    )
                else:
                    for record_item in provideritem.record.recorditem:
                        if record_item.fileid not in fileids:
                            fileids.add(record_item.fileid)
                            providers[provider].record.recorditem.append(
                                record_item
                            )
                            providers[provider].no_of_records_found += 1
                            providers[provider].no_of_records_returned += 1
        response = self.api.factory.create('QueryResponse')
        response.provideritem = providers.values()
        return response
    
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
            downloader = download.Downloader()
            threading.Thread(target=downloader.reactor.run).start()
            res = Results(
                lambda _: downloader.reactor.stop(), 1,
                lambda mp: self.link(query_response, mp)
            )
        else:
            res = Results(
                lambda _: None, 1, lambda mp: self.link(qr, mp)
            )
        if path is None:
            path = os.path.join(tempfile.mkdtemp(), '{file}')
        self.download_all(
            self.api.service.GetData(
                self.make_getdatarequest(query_response, methods)
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
                item = _Str(map_[record_item.fileid]['path'])
                item.meta = record_item
                ret.append(item)
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
                    try:
                        self.download(
                            dresponse.method.methodtype[0],
                            dataitem.url,
                            dw,
                            res.require(map(str, dataitem.fileiditem.fileid)),
                            path,
                            qr[dataitem.fileiditem.fileid[0]]
                        )
                    except ValueError:
                        # TODO: Log
                        continue
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
                    try:
                        info = self.missing_information(
                            info, dresponse.info
                        )
                    except ValueError:
                        # TODO: Log.
                        continue
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
        raise ValueError
    
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
        raise ValueError


class InteractiveAPI(API):
    def multiple_choices(self, choices, response):
        while True:
            for n, elem in enumerate(choices):
                print "(%d) %s" % (n + 1, elem)
            choice = raw_input("Method number: ")
            try:
                return [choices[int(choice) - 1]]
            except ValueError, IndexError:
                continue
            except KeyboardInterrupt:
                raise ValueError
        
    def missing_information(self, info, field):
        return raw_input(field + ': ')


if __name__ == '__main__':
    import sunpy
    api = API()
    
    qr = api.query_legacy(
        datetime(2010, 1, 1), datetime(2010, 1, 1, 1),
        instrument='eit'
    )
    
    res = api.get(qr).wait()
    sunpy.Map(res[0]).plot()
