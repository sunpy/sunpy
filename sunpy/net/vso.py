# -*- coding: utf-8 -*-
# Author: Florian Mayer <florian.mayer@bitsrc.org>

import re
import os
import sys
import tempfile
import threading

from datetime import datetime, timedelta
from functools import partial
from collections import defaultdict

from suds import client

from sunpy.net import download
from sunpy.util.util import anytim

DEFAULT_URL = 'http://docs.virtualsolar.org/WSDL/VSOi_rpc_literal.wsdl'
DEFAULT_PORT = 'nsoVSOi'
TIMEFORMAT = '%Y%m%d%H%M%S'
RANGE = re.compile(r'(\d+)(\s*-\s*(\d+))?(\s*([a-zA-Z]+))?')


# TODO: Name
class NoData(Exception):
    pass


class _Str(str):
    pass


class _Attr(object):
    def __and__(self, other):
        if isinstance(other, _AttrOr):
            return _AttrOr([elem & self for elem in other.attrs])
        if isinstance(other, self.__class__):
            # A record cannot match two different values
            # for the same attribute.
            # TODO: Error?
            return NotImplemented
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
        for elem in self.attrs:
            if isinstance(other, elem.__class__):
                # A record cannot match two different values
                # for the same attribute.
                # TODO: Error?
                return NotImplemented
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
        value = api.factory.create('QueryRequestBlock')
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
        value = api.factory.create('QueryRequestBlock')
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
            return (min_ / k, max_ / k)
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
        self.errors = []
    
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
    
    def add_error(self, exception):
        self.errors.append(exception)


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


def _parse_waverange(string):
    min_, max_, unit = RANGE.match(string)[::2]
    return {
        'wave_wavemin': min_,
        'wave_wavemax': min_ if max_ is None else max_,
        'wave_waveunit': 'Angstrom' if unit is None else unit,
    }


def _parse_date(string):
    start, end = string.split(' - ')
    return {'time_start': start.strip(), 'time_end': end.strip()}


def iter_records(response):
    for prov_item in response.provideritem:
        if not hasattr(prov_item, 'record') or not prov_item.record:
            continue
        for record_item in prov_item.record.recorditem:
            yield record_item


def iter_errors(response):
    for prov_item in response.provideritem:
        if not hasattr(prov_item, 'record') or not prov_item.record:
            yield prov_item


class QueryResponse(list):
    def __init__(self, lst, queryresult=None):
        super(QueryResponse, self).__init__(lst)
        self.queryresult = queryresult
    
    @classmethod
    def create(cls, queryresult):
        return cls(iter_records(queryresult), queryresult)
    
    def total_size(self):
        # Warn about -1 values?
        return sum(record.size for record in self if record.size > 0)
    
    def no_records(self):
        return sum(1 for _ in self)
    
    def time_range(self):
        return (
            datetime.strptime(
                min(record.time.start for record in self), TIMEFORMAT),
            datetime.strptime(
                max(record.time.end for record in self), TIMEFORMAT)
        )


class DownloadFailed(Exception):
    pass

class MissingInformation(Exception):
    pass

class UnknownMethod(Exception):
    pass

class MultipleChoices(Exception):
    pass

class UnknownVersion(Exception):
    pass

class UnknownStatus(Exception):
    pass

class VSOClient(object):
    method_order = [
        'URL-TAR_GZ', 'URL-ZIP', 'URL-TAR', 'URL-FILE', 'URL-packaged'
    ]
    def __init__(self, api=None):
        if api is None:
            api = client.Client(DEFAULT_URL)
            api.set_options(port=DEFAULT_PORT)
        self.api = api
    
    def make(self, type_, **kwargs):
        obj = self.api.factory.create(type_)
        for k, v in kwargs.iteritems():
            split = k.split('__')
            tip = split[-1]
            rest = split[:-1]
            
            item = obj
            for elem in rest:
                item = item[elem]
            
            if isinstance(v, dict):
                # Do not throw away type information for dicts.
                for k, v in v.iteritems():
                    item[tip][k] = v
            else:
                item[tip] = v
        return obj
    
    def query(self, query):
        return QueryResponse.create(self.merge(
            self.api.service.Query(self.make('QueryRequest', block=block))
            for block in query.create(self.api)
        ))
    
    def merge(self, queryresponses):      
        fileids = set()
        providers = {}
        
        for queryresponse in queryresponses:
            for provideritem in queryresponse.provideritem:
                provider = provideritem.provider
                if not hasattr(provideritem.record, 'recorditem'):
                    continue
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
        return self.make('QueryResponse', provideritem=providers.values())
    
    def query_legacy(self, tstart=None, tend=None, **kwargs):
        sdk = lambda key: lambda value: {key: value}
        ALIASES = {
            'wave_min': sdk('wave_wavemin'),
            'wave_max': sdk('wave_wavemax'),
            'wave_type': sdk('wave_wavetype'),
            'wave_unit': sdk('wave_waveunit'),
            'wave': _parse_waverange,
            'inst': sdk('instrument'),
            'telescope': sdk('instrument'),
            'spacecraft': sdk('source'),
            'observatory': sdk('source'),
            'start_date': sdk('time_start'),
            'end_date': sdk('time_end'),
            'start': sdk('time_start'),
            'end': sdk('time_end'),
            'date': _parse_date,
        }
        kwargs.update({'time_start': tstart, 'time_end': tend})
        
        queryreq = self.api.factory.create('QueryRequest')
        for key, value in kwargs.iteritems():
            if key.startswith('time'):
                value = anytim(value).strftime(TIMEFORMAT)
            for k, v in ALIASES.get(key, sdk(key))(value).iteritems():
                attr = k.split('_')
                lst = attr[-1]
                rest = attr[:-1]
                
                item = queryreq.block
                for elem in rest:
                    try:
                        item = item[elem]
                    except KeyError:
                        raise ValueError("Unexpected argument %s." % key)
                if lst not in item:
                    raise ValueError("Unexpected argument %s." % key)
                if item[lst]:
                    raise ValueError("Got multiple values for %s." % k)
                item[lst] = v
        return QueryResponse.create(self.api.service.Query(queryreq))
    
    def latest(self):
        return self.query_legacy(
            datetime.utcnow()  - timedelta(7),
            datetime.utcnow(),
            time_near=datetime.utcnow()
        )
    
    def get(self, query_response, path=None, methods=['URL-FILE'], downloader=None):
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
        fileids = VSOClient.by_fileid(query_response)
        if not fileids:
            res.poke()
            return res
        self.download_all(
            self.api.service.GetData(
                self.make_getdatarequest(query_response, methods)
                ),
            methods, downloader, path,
            VSOClient.by_fileid(query_response), res
        )
        res.poke()
        return res
    
    @staticmethod
    def link(query_response, map_):
        if not map_:
            return []
        ret = []
        
        for record_item in query_response:
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
        
        return self.make(
            'VSOGetDataRequest',
            request__method__methodtype=methods, 
            request__info=info,
            request__datacontainer__datarequestitem=[
                self.make('DataRequestItem', provider=k, fileiditem__fileid=[v])
                for k, v in map_.iteritems()
            ]
        )
    
    def download_all(self, response, methods, dw, path, qr, res, info=None):
        GET_VERSION = [
            ('0.8', (5, 8)),
            ('0.7', (1, 4)),
            ('0.6', (0, 3)),
        ]
        for dresponse in response.getdataresponseitem:
            for version, (from_, to) in GET_VERSION:
                if dresponse.version >= version:
                    break
            else:
                res.add_error(UnknownVersion(dresponse))
                continue
            
            code = (
                dresponse.status[from_:to]
                if hasattr(dresponse, 'status') else '200'
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
                    except NoData:
                        res.add_error(DownloadFailed(dresponse))
                        continue
            elif code == '300' or code == '412' or code == '405':
                if code == '300':
                    try:
                        methods = self.multiple_choices(
                            dresponse.method.methodtype, dresponse
                        )
                    except NoData:
                        res.add_error(MultipleChoices(dresponse))
                        continue
                elif code == '412':
                    try:
                        info = self.missing_information(
                            info, dresponse.info
                        )
                    except NoData:
                        res.add_error(MissingInformation(dresponse))
                        continue
                elif code == '405':
                    try:
                        methods = self.unknown_method(dresponse)
                    except NoData:
                        res.add_error(UnknownMethod(dresponse))
                        continue
                
                files = []
                for dataitem in dresponse.getdataitem.dataitem:
                    files.extend(dataitem.fileiditem.fileid)
                
                request = self.create_getdatarequest(
                    {dresponse.provider: files}, methods, info
                )
                
                self.download_all(
                    self.api.service.GetData(request), methods, dw, path,
                    qr, res, info
                )
            else:
                res.add_error(UnknownStatus(dresponse))
    
    def download(self, method, url, dw, callback, *args):
        if method.startswith('URL'):
            dw.reactor.call_sync(
                partial(dw.download, url, partial(mk_filename, *args),
                        callback)
            )
        raise NoData
    
    @staticmethod
    def by_provider(response):
        map_ = defaultdict(list)
        for record in response:
            map_[record.provider].append(record)
        return map_
    
    @staticmethod
    def by_fileid(response):
        return dict(
            (record.fileid, record) for record in response
        )
    
    def multiple_choices(self, choices, response):
        for elem in self.method_order:
            if elem in choices:
                return [elem]
        raise NoData
    
    def missing_information(self, info, field):
        raise NoData
    
    def unknown_method(self, response):
        raise NoData


class InteractiveVSOClient(VSOClient):
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
                raise NoData
        
    def missing_information(self, info, field):
        return raw_input(field + ': ')
    
    def search(self, tstart=None, tend=None, **kwargs):
        if isinstance(tstart, _Attr):
            return self.query(tstart)
        else:
            return self.query_legacy(tstart, tend, **kwargs)


g_client = None
def search(tstart=None, tend=None, **kwargs):
    """ When passed a single _Attr object, perform new-style query;
    otherwise, perform legacy query.
    """
    global g_client
    if g_client is None:
        g_client = InteractiveVSOClient()
    return g_client.search(tstart, tend, **kwargs)


def get(query_response, path=None, methods=['URL-FILE'], downloader=None):
    """
    Download data specified in the query_response.
    
    Parameters
    ----------
    query_response : sunpy.net.vso.QueryResponse
        QueryResponse containing the items to be downloaded.
    path : str
        Specify where the data is to be downloaded. Can refer to arbitrary
        fields of the QueryResponseItem (instrument, source, time, ...) via
        string formatting, moreover the file-name of the file downloaded can
        be refered to as file, e.g.
        "{source}/{instrument}/{time.start}/{file}".
    methods : {list of str}
        Methods acceptable to user.
    downloader : sunpy.net.downloader.Downloader
        Downloader used to download the data.
    
    """
    global g_client
    if g_client is None:
        g_client = InteractiveVSOClient()
    return g_client.get(query_response, path, methods, downloader)


# Add latest?
if __name__ == '__main__':
    import sunpy
    qr = search(
        datetime(2010, 1, 1), datetime(2010, 1, 1, 1),
        instrument='eit'
    )
    
    res = get(qr).wait()
    
    print res[0]
    sunpy.Map(res[0]).plot()
