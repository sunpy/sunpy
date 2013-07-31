from sunpy.time import parse_time
from sunpy.net.vso import attrs as vso_attrs
from sunpy.net.attr import AttrWalker, Attr, ValueAttr, AttrAnd, AttrOr
from sunpy.database.tables import DatabaseEntry, Tag as TableTag

__all__ = ['Starred', 'Tag', 'walker']


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
        return '<%s(%s)>' % (self.__class__.__name__, self.value)


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
        return '<Tag(%r, %r)>' % (self.tagname, self.inverted)


class Path(vso_attrs._VSOSimpleAttr):
    pass


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
        download_time .inverted = True
        return download_time

    def collides(self, other):  # pragma: no cover
        return isinstance(other, self.__class__)

    def __repr__(self):
        return '<%sDownloadTime(%r, %r)>' % (
            '~' if self.inverted else '', self.start, self.end)


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
        elif typ == 'download time':
            start, end, inverted = value
            if inverted:
                query = query.filter(~DatabaseEntry.download_time.between(start, end))
            else:
                query = query.filter(DatabaseEntry.download_time.between(start, end))
        else:
            query = query.filter_by(**{typ: value})
    return query.all()


@walker.add_converter(Tag)
def _convert(attr):
    return ValueAttr({('tag', attr.inverted): attr.tagname})


@walker.add_converter(Starred)
def _convert(attr):
    return ValueAttr({('starred', ): attr.value})


@walker.add_converter(vso_attrs._VSOSimpleAttr)
def _convert(attr):
    return ValueAttr({(attr.__class__.__name__.lower(),): attr.value})


@walker.add_converter(DownloadTime)
def _convert(attr):
    return ValueAttr({
        ('download time', ): (attr.start, attr.end, attr.inverted)})
