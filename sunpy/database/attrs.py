from sqlalchemy import and_, not_, or_

from sunpy.database.tables import DatabaseEntry
from sunpy.database.tables import FitsHeaderEntry as TableFitsHeaderEntry
from sunpy.database.tables import Tag as TableTag
from sunpy.net import _attrs as core_attrs
from sunpy.net.attr import Attr, AttrAnd, AttrOr, AttrWalker, Range, SimpleAttr, ValueAttr
from sunpy.time import parse_time

__all__ = [
    'Starred', 'Tag', 'Path', 'DownloadTime', 'FitsHeaderEntry', 'walker']

# This frozenset has been hardcoded to denote VSO attributes that are
# currently supported, on derdon's request.
SUPPORTED_SIMPLE_VSO_ATTRS = frozenset(['source', 'provider', 'physobs',
                                        'instrument'])
SUPPORTED_NONVSO_ATTRS = frozenset(['starred'])


class _BooleanAttr:
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

    def __nonzero__(self):  # py 2.x
        return self.value

    def __bool__(self):  # py 3.x
        return self.value

    def __invert__(self):
        attr = self.make()
        attr.value = not self.value
        return attr

    def __eq__(self, other):
        return isinstance(other, self.make) and self.value == other.value

    def __hash__(self):
        return super().__hash__()

    def collides(self, other):  # pragma: no cover
        return False

    def __repr__(self):
        return '<{}{}()>'.format(
            '~' if not self.value else '', self.__class__.__name__)


class Starred(_BooleanAttr, Attr):
    type_name = "starred"

    def __init__(self):
        super().__init__(True, self.__class__)


class Tag(Attr):
    type_name = "tag"

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
        return '<{}Tag({!r})>'.format(
            '~' if self.inverted else '', self.tagname)


class Path(Attr):
    type_name = 'path'

    def __init__(self, value, inverted=False):
        self.value = value
        self.inverted = inverted

    def __invert__(self):
        return self.__class__(self.value, True)

    def collides(self, other):  # pragma: no cover
        return isinstance(other, self.__class__)

    def __repr__(self):
        return '<{}Path({!r})>'.format(
            '~' if self.inverted else '', self.value)


# TODO: support excluding ranges as soon as
# attr.Range.__xor__ is fixed / renamed
class DownloadTime(Range):
    type_name = 'download time'

    def __init__(self, start, end):
        self.start = parse_time(start).datetime
        self.end = parse_time(end).datetime
        self.inverted = False
        super().__init__(start, end)

    def __invert__(self):
        download_time = self.__class__(self.start, self.end)
        download_time.inverted = True
        return download_time

    def collides(self, other):  # pragma: no cover
        return isinstance(other, self.__class__)

    def __repr__(self):
        return '<{}DownloadTime({!r}, {!r})>'.format(
            '~' if self.inverted else '', self.start, self.end)


class FitsHeaderEntry(Attr):
    type_name = "fitsheaderentry"

    def __init__(self, key, value, inverted=False):
        self.key = key
        self.value = value
        self.inverted = inverted

    def __invert__(self):
        return self.__class__(self.key, self.value, True)

    def collides(self, other):  # pragma: no cover
        return False

    def __repr__(self):
        return '<{}FitsHeaderEntry({!r}, {!r})>'.format(
            '~' if self.inverted else '', self.key, self.value)


walker = AttrWalker()


@walker.add_creator(AttrOr)
def _create(wlk, root, session):
    entries = [set(wlk.create(attr, session)) for attr in root.attrs]
    return list(set.union(*entries))


@walker.add_creator(AttrAnd)
def _create(wlk, root, session):
    entries = [set(wlk.create(attr, session)) for attr in root.attrs]
    return list(set.intersection(*entries))


def _inverter_helper(query, inverted):
    return not_(query) if inverted else query


@walker.add_creator(ValueAttr)
def _create(wlk, root, session):
    query = session.query(DatabaseEntry)
    for key, value in root.attrs.items():
        typ = key[0].lower()
        if typ == Tag.type_name:
            criterion = TableTag.name.in_([value])
            base_query = DatabaseEntry.tags.any(criterion)
            # `key[1]` is here the `inverted` attribute of the tag. That means
            # that if it is True, the given tag must not be included in the
            # resulting entries.
            query = query.filter(_inverter_helper(base_query, key[1]))

        elif typ == FitsHeaderEntry.type_name:
            key, val, inverted = value
            key_criterion = TableFitsHeaderEntry.key == key
            value_criterion = TableFitsHeaderEntry.value == val
            base_query = and_(
                DatabaseEntry.fits_header_entries.any(key_criterion),
                DatabaseEntry.fits_header_entries.any(value_criterion))
            query = query.filter(_inverter_helper(base_query, inverted))

        elif typ == DownloadTime.type_name:
            start, end, inverted = value
            base_query = DatabaseEntry.download_time.between(start, end)
            query = query.filter(_inverter_helper(base_query, inverted))

        elif typ == Path.type_name:
            path, inverted = value
            base_query = _inverter_helper(DatabaseEntry.path == path, inverted)
            if inverted:
                base_query = or_(base_query, DatabaseEntry.path == None)  # NOQA
            query = query.filter(base_query)

        elif typ == core_attrs.Wavelength.type_name:
            wavemin, wavemax, waveunit = value
            query = query.filter(and_(
                DatabaseEntry.wavemin >= wavemin,
                DatabaseEntry.wavemax <= wavemax))

        elif typ == core_attrs.Time.type_name:
            start, end, _ = value
            query = query.filter(and_(
                DatabaseEntry.observation_time_start < end,
                DatabaseEntry.observation_time_end > start))

        elif typ in (SUPPORTED_SIMPLE_VSO_ATTRS | SUPPORTED_NONVSO_ATTRS):
            query = query.filter_by(**{typ: value})

        else:
            raise NotImplementedError(
                f"The attribute {typ!r} is not yet supported to query a database.")

    return query.all()


@walker.add_converter(Tag)
def _convert(attr):
    return ValueAttr({(Tag.type_name, attr.inverted): attr.tagname})


@walker.add_converter(Starred)
def _convert(attr):
    return ValueAttr({(Starred.type_name, ): attr.value})


@walker.add_converter(Path)
def _convert(attr):
    return ValueAttr({(Path.type_name, ): (attr.value, attr.inverted)})


@walker.add_converter(DownloadTime)
def _convert(attr):
    return ValueAttr({
        (DownloadTime.type_name, ): (attr.start, attr.end, attr.inverted)})


@walker.add_converter(FitsHeaderEntry)
def _convert(attr):
    return ValueAttr(
        {(FitsHeaderEntry.type_name, ): (attr.key, attr.value, attr.inverted)})


@walker.add_converter(SimpleAttr)
def _convert(attr):
    return ValueAttr({(attr.type_name, ): attr.value})


@walker.add_converter(core_attrs.Wavelength)
def _convert(attr):
    return ValueAttr({(attr.type_name, ): (attr.min.value, attr.max.value, str(attr.unit))})


@walker.add_converter(core_attrs.Time)
def _convert(attr):
    near = None if not attr.near else attr.near.datetime
    return ValueAttr({(attr.type_name, ): (attr.start.datetime, attr.end.datetime, near)})
