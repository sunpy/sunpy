# Author: Simon Liedtke <liedtke.simon@googlemail.com>
#
# This module was developed with funding provided by
# the Google Summer of Code (2013).

from __future__ import absolute_import

from sqlalchemy import or_, and_, not_

from sunpy.time import parse_time
from sunpy.net.vso import attrs as vso_attrs
from sunpy.net.attr import AttrWalker, Attr, ValueAttr, AttrAnd, AttrOr
from sunpy.database.tables import DatabaseEntry, Tag as TableTag,\
    FitsHeaderEntry as TableFitsHeaderEntry

__all__ = [
    'Starred', 'Tag', 'Path', 'DownloadTime', 'FitsHeaderEntry', 'walker']

# This frozenset has been hardcoded to denote VSO attributes that are currently supported, on derdon's request.
SUPPORTED_SIMPLE_VSO_ATTRS = frozenset(['source', 'provider', 'physobs', 'instrument'])
SUPPORTED_NONVSO_ATTRS = frozenset(['starred'])

class _BooleanAttr(object):
    def __init__(self, value, make):
        self.value = bool(value)
        self.make = make

    def __and__(self, other):
        if not isinstance(other, self.make):
            return AttrAnd([self, other])
        attr = self.make()
        attr.value = self.value and other.value
        return attr

    def __or__(self, other):
        if not isinstance(other, self.make):
            return AttrOr([self, other])
        attr = self.make()
        attr.value = self.value or other.value
        return attr

    def __nonzero__(self):
        return self.value

    def __invert__(self):
        attr = self.make()
        attr.value = not self.value
        return attr

    def __eq__(self, other):
        return isinstance(other, self.make) and self.value == other.value

    def collides(self, other):  # pragma: no cover
        return False

    def __repr__(self):
        return '<{0}{1}()>'.format(
            '~' if not self.value else '', self.__class__.__name__)


class Starred(_BooleanAttr, Attr):
    def __init__(self):
        _BooleanAttr.__init__(self, True, self.__class__)


class Tag(Attr):
    def __init__(self, tagname):
        self.tagname = tagname
        self.inverted = False

    def __invert__(self):
        tag = self.__class__(self.tagname)
        tag.inverted = True
        return tag

    def collides(self, other):  # pragma: no cover
        return False

    def __repr__(self):
        return '<{0}Tag({1!r})>'.format(
            '~' if self.inverted else '', self.tagname)


class Path(Attr):
    def __init__(self, value):
        self.value = value
        self.inverted = False

    def __invert__(self):
        path = self.__class__(self.value)
        path.inverted = True
        return path

    def collides(self, other):  # pragma: no cover
        return isinstance(other, self.__class__)

    def __repr__(self):
        return '<{0}Path({1!r})>'.format(
            '~' if self.inverted else '', self.value)


# TODO: support excluding ranges as soon as
# vso_attrs._Range.__xor__ is fixed / renamed
class DownloadTime(Attr, vso_attrs._Range):
    def __init__(self, start, end):
        self.start = parse_time(start)
        self.end = parse_time(end)
        self.inverted = False
        vso_attrs._Range.__init__(self, start, end, self.__class__)

    def __invert__(self):
        download_time = self.__class__(self.start, self.end)
        download_time.inverted = True
        return download_time

    def collides(self, other):  # pragma: no cover
        return isinstance(other, self.__class__)

    def __repr__(self):
        return '<{0}DownloadTime({1!r}, {2!r})>'.format(
            '~' if self.inverted else '', self.start, self.end)


class FitsHeaderEntry(Attr):
    def __init__(self, key, value):
        self.key = key
        self.value = value
        self.inverted = False

    def __invert__(self):
        entry = self.__class__(self.key, self.value)
        entry.inverted = True
        return entry

    def collides(self, other):  # pragma: no cover
        return False

    def __repr__(self):
        return '<{0}FitsHeaderEntry({1!r}, {2!r})>'.format(
            '~' if self.inverted else '', self.key, self.value)


walker = AttrWalker()


@walker.add_creator(AttrOr)
def _create(wlk, root, session):
    entries = set()
    for attr in root.attrs:
        # make sure to add only new entries to the set to avoid duplicates
        entries.update(set(wlk.create(attr, session)) - entries)
    return list(entries)


@walker.add_creator(AttrAnd)
def _create(wlk, root, session):
    entries = [set(wlk.create(attr, session)) for attr in root.attrs]
    return list(set.intersection(*entries))


@walker.add_creator(ValueAttr)
def _create(wlk, root, session):
    query = session.query(DatabaseEntry)
    for key, value in root.attrs.iteritems():
        typ = key[0]
        if typ == 'tag':
            criterion = TableTag.name.in_([value])
            # `key[1]` is here the `inverted` attribute of the tag. That means
            # that if it is True, the given tag must not be included in the
            # resulting entries.
            if key[1]:
                query = query.filter(~DatabaseEntry.tags.any(criterion))
            else:
                query = query.filter(DatabaseEntry.tags.any(criterion))
        elif typ == 'fitsheaderentry':
            key, val, inverted = value
            key_criterion = TableFitsHeaderEntry.key == key
            value_criterion = TableFitsHeaderEntry.value == val
            if inverted:
                query = query.filter(not_(and_(
                    DatabaseEntry.fits_header_entries.any(key_criterion),
                    DatabaseEntry.fits_header_entries.any(value_criterion))))
            else:
                query = query.filter(and_(
                    DatabaseEntry.fits_header_entries.any(key_criterion),
                    DatabaseEntry.fits_header_entries.any(value_criterion)))
        elif typ == 'download time':
            start, end, inverted = value
            if inverted:
                query = query.filter(
                    ~DatabaseEntry.download_time.between(start, end))
            else:
                query = query.filter(
                    DatabaseEntry.download_time.between(start, end))
        elif typ == 'path':
            path, inverted = value
            if inverted:
                # pylint: disable=E711
                query = query.filter(or_(
                    DatabaseEntry.path != path, DatabaseEntry.path == None))
            else:
                query = query.filter(DatabaseEntry.path == path)
        elif typ == 'wave':
            wavemin, wavemax, waveunit = value
            query = query.filter(and_(
                DatabaseEntry.wavemin >= wavemin,
                DatabaseEntry.wavemax <= wavemax))
        elif typ == 'time':
            start, end, near = value
            query = query.filter(and_(
                DatabaseEntry.observation_time_start < end,
                DatabaseEntry.observation_time_end > start))
        else:
            if typ.lower() not in SUPPORTED_SIMPLE_VSO_ATTRS.union(SUPPORTED_NONVSO_ATTRS):
                raise NotImplementedError("The attribute {0!r} is not yet supported to query a database.".format(typ))
            query = query.filter_by(**{typ: value})
    return query.all()


@walker.add_converter(Tag)
def _convert(attr):
    return ValueAttr({('tag', attr.inverted): attr.tagname})


@walker.add_converter(Starred)
def _convert(attr):
    return ValueAttr({('starred', ): attr.value})


@walker.add_converter(Path)
def _convert(attr):
    return ValueAttr({('path', ): (attr.value, attr.inverted)})


@walker.add_converter(DownloadTime)
def _convert(attr):
    return ValueAttr({
        ('download time', ): (attr.start, attr.end, attr.inverted)})


@walker.add_converter(FitsHeaderEntry)
def _convert(attr):
    return ValueAttr(
        {('fitsheaderentry', ): (attr.key, attr.value, attr.inverted)})


@walker.add_converter(vso_attrs._VSOSimpleAttr)
def _convert(attr):
    return ValueAttr({(attr.__class__.__name__.lower(), ): attr.value})


@walker.add_converter(vso_attrs.Wave)
def _convert(attr):
    return ValueAttr({('wave', ): (attr.min.value, attr.max.value, str(attr.unit))})


@walker.add_converter(vso_attrs.Time)
def _convert(attr):
    return ValueAttr({('time', ): (attr.start, attr.end, attr.near)})
