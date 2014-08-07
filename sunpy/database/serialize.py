from operator import attrgetter
from datetime import datetime
import json
from astropy import units as u

from sunpy.net import vso
from sunpy.database import attrs as db_attrs
from sunpy.net.attr import Attr, AttrOr, AttrAnd


__all__ = ['dump_query', 'load_query']


class QueryEncoder(json.JSONEncoder):
    def default(self, o):
        if isinstance(o, datetime):
            return o.strftime("%Y-%m-%d %H:%M:%S")
        elif isinstance(o, Attr):
            if isinstance(o, (AttrAnd, AttrOr)):
                # sort by dictionary keys to be order-invariant
                values = sorted(o.attrs, key=attrgetter('__class__.__name__'))
            elif isinstance(o, vso.attrs.Wave):
                values = o.min.value, o.max.value, str(o.unit)
            elif isinstance(o, vso.attrs.Time):
                values = o.start, o.end, o.near
            elif isinstance(o, (vso.attrs._VSOSimpleAttr, db_attrs.Starred)):
                values = o.value
            elif isinstance(o, db_attrs.Tag):
                values = o.tagname, o.inverted,
            elif isinstance(o, db_attrs.Path):
                values = o.value, o.inverted,
            elif isinstance(o, db_attrs.DownloadTime):
                values = o.start, o.end, o.inverted,
            elif isinstance(o, db_attrs.FitsHeaderEntry):
                values = o.key, o.value, o.inverted,
            else:
                assert False
            return {o.__class__.__name__: values}
        return json.JSONEncoder.default(self, o)


def query_decode(json_object):
    simple_attrs = [
        'Provider', 'Source', 'Instrument', 'Physobs', 'Pixels', 'Level',
        'Resolution', 'Detector', 'Filter', 'Sample', 'Quicklook', 'PScale']
    for key in simple_attrs:
        if key in json_object:
            Attr = getattr(vso.attrs, key)
            return Attr(json_object[key])
    if 'Wave' in json_object:
        Attr = getattr(vso.attrs, 'Wave')
        wavemin, wavemax, unit = json_object['Wave']
        return Attr(wavemin * u.Unit(unit), wavemax * u.Unit(unit))
    if 'Time' in json_object:
        Attr = getattr(vso.attrs, 'Time')
        return Attr(*json_object['Time'])
    for key in ['Tag', 'Path', 'DownloadTime', 'FitsHeaderEntry']:
        if key in json_object:
            Attr = getattr(db_attrs, key)
            values, inverted = json_object[key][:-1], json_object[key][-1]
            if inverted:
                return ~Attr(*values)
            return Attr(*values)
    if 'Starred' in json_object:
        if json_object['Starred']:
            return ~db_attrs.Starred()
        return db_attrs.Starred()
    for key in ['AttrOr', 'AttrAnd']:
        if key in json_object:
            Attr = getattr(vso.attrs, key)
            return Attr(json_object[key])
    return json_object


def dump_query(query):  # pragma: no cover
    return json.dumps(query, cls=QueryEncoder)


def load_query(dump):  # pragma: no cover
    return json.loads(dump, object_hook=query_decode)
